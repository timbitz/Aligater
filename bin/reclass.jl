#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 3/16/2015
   e-mail: tim.sterne.weiler@utoronto.ca
 
   julia 0.4-dev
=# 

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ## INITIALIZATION  ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

global head = "[aligater reclass]:"
println(STDERR, "$head Loading Packages..")

using ArgParse
using StatsBase
using Match

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--geneid", "-g"
      help = "flag to collapse geneids"
      action = :store_true
    "--biotype", "-b"
      help = "flag to collapse biotypes"
      action = :store_true
    "--uniq", "-u"
      help = "filter for unique chimeras by left/right lengths and 10bp unique match of lig site"
      action = :store_true
    "--load"
      help = "previously stored background file for 'dclass' in HDF5 format"
      arg_type = ASCIIString
    "--save"
      help = "store to file the resulting 'dclass' in jlz format as calculated on the input data"
      arg_type = ASCIIString 
    "--barcode"
      help = "random barcode associated with each read in jlz format (optional)"
      arg_type = ASCIIString
    "--size"
      help = "nmer size to use for index-based reclassification"
      arg_type = Int
      default = 32
  end
  return parse_args(s)
end

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ###  FUNCTIONS ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

include("dictext.jl") #--> savedict, loaddict, dush!, dinc!, dnorm!

function isUniqueJunc!{K <: AbstractString}( used::Dict{K,Bool}, seq, genes, dbar::Dict{K,K}, name )
  barcode = get(dbar, name, "")
  m = match(r"([AGCTUN]{5}_[AGCTUN]{5})", seq)
  cap = m.captures[1]
  m = match(r"([AGCTUN]+)_([AGCTUN]+)", seq)
  lLen,rLen = length(m.captures[1]), length(m.captures[2])
  geneid = sort( genes )
  key = join(geneid, ":") * "_$cap:$lLen:" * barcode
  retbool = false
  if !( haskey( used, key ) )
      used[key] = true
      retbool = true
  end
  retbool
end #--> Bool


# This function can both be used to set, and also get the class associated with two chimeric
# sequences based on window size homology.  First pass should set the database then run a
# second pass to retrieve the final results
function setClass!( dict::Dict{ASCIIString,Char}, class::Char, seqs; size=32, stepsize=floor(size/3), set=false )
  @assert( length(seqs) >= 2 )
  classHeirPair( a::Char, b::Char ) = a < b ? b : a  #--> Char
  classHeirArray( arr::Array{Char,1} ) = max( arr... ) #--> Char
#  retarray = Char[]
  retval = class
  lenA,lenB = map(length, seqs)
  # pad a sequence if it is shorter than our target window size
  seqA = seqs[1] * repeat(".", max( size - lenA, 0 ))
  seqB = seqs[2] * repeat(".", max( size - lenB, 0 ))
  # reset the lengths after padding
  lenA,lenB = length(seqA),length(seqB)
  it=max(1, Int(stepsize)) # calculate the step size, as half the window size
  # iterate through the cartesian product of windows in seqA and seqB by it step size
  for i in 1:it:(lenA-size)+1, j in 1:it:(lenB-size)+1
    winA = seqA[i:(i+size-1)] # access substrings 
    winB = seqB[j:(j+size-1)]
    keyA,keyB = winA < winB ? (winA,winB) : (winB,winA)
    jointkey = keyA * "_" * keyB
    # use non-exported ht_keyindex function to avoid over indexing hash
    indx = Base.ht_keyindex( dict, jointkey )
    if indx <= 0
      if set
        dict[jointkey] = class #set val if none exists
        #length(retarray) <= 0 && push!(retarray, class)
      end
    else
      curVal = classHeirPair( dict.vals[indx], class )
      dict.vals[indx] = curVal
      #push!(retarray, curVal)
      retval = classHeirPair( retval, curVal )
    end
  end
#  @assert( length(retarray) > 0 )
#  length(retarray) > 1 ? classHeirArray( retarray ) : retarray[1]
  retval
end #--> Char

function reducegeneid( geneid, biotype, repeatname, repeatclass )
  retval = @match repeatname begin
             r"SSU-rRNA" => "18S_rRNA"
             r"LSU-rRNA" => "28S_rRNA"
             _, if ismatch(r"RNA", repeatclass) end => repeatname
             _ => geneid
           end
   (retval != geneid) && return retval

  retval = @match geneid begin
             r"RNA{0,1}5S|C15orf52" => "5S_rRNA"
             r"RNA{0,1}5-8" => "5.8S_rRNA"
             r"7SL" => "7SLRNA_srpRNA"
             r"7SK" => "7SK_RNA"
             _ => geneid
          end
  (retval != geneid) && return retval

  m = match(r"RNU(\d+)(?=(ATAC|-))", geneid)
  if m != nothing
    retval = @match m.captures[2] begin
               "ATAC" => "RNU"*m.captures[1]*m.captures[2]
               "-" => "U"*m.captures[1]*"_snRNA"
               _ => geneid
             end
  end
  (retval != geneid) && return retval

  m = match(r"RNA{0,1}(\d+)S", geneid)
  if m != nothing
    retval = m.captures[1] * "S_rRNA"
  end
  retval
end #--> ASCIIString

function reducebiotype( geneid, biotype, repeatname, repeatclass )
  biotype = replace(biotype, r"Mt-", "")
  ismatch(r"(srp|sn|r)RNA", biotype) && return(biotype)
  ismatch(r"tRNA", repeatclass) && return("tRNA")
  ismatch(r"SNOR|SCAR|^ACA", geneid) && return("snoRNA")
  ismatch(r"RMRP", geneid) && return("RNase_MRP")
  ismatch(r"RPPH", geneid) && return("RNase_P")
  ismatch(r"VT.*RNA", geneid) && return("vtRNA")
  ismatch(r"7SK_RNA", repeatname) && return(repeatname)
  ismatch(r"U(3|8|9|\d{2,+})|snoRNA", repeatname) && return("snoRNA")   
  ismatch(r"snRNA", repeatname) && return("snRNA")
  if repeatclass != "NA"
    return("$biotype\_$repeatclass")
  else
    return(biotype == "NA" ? geneid : biotype)
  end
end #--> ASCIIString

function reducef( func, str )
  const geneInd = 4 # gene-id index
  const biotInd = 7 # biotype index
  const repnInd = 8 # repeat name index
  const repcInd = 9 # repeat class index
  
  function safesplit( arr, ind; char=':' )
    @assert( length(arr) >= ind )
    split( arr[ind], char )
  end

  gspl = safesplit( str, geneInd )
  bspl = safesplit( str, biotInd )
  nspl = safesplit( str, repnInd )
  cspl = safesplit( str, repcInd )

  lenLim = min( map(length, (gspl, bspl, nspl, cspl))... )
  res = ASCIIString[]
  for i in 1:lenLim
    push!(res, func(gspl[i], bspl[i], nspl[i], cspl[i]))
  end
  res
end #--> ASCIIString[]

function reclassReduceAndPrint( dclass::Dict, dstore::Dict, pargs, seqInd )
  for class in keys(dstore), s in dstore[class]
    seqs = masksplit(s[seqInd], '_')
    if length(seqs) > 1 # try to reclassify
      s[1] = convert(ASCIIString, string( setClass!( dclass, s[1][1], seqs, stepsize=1 ) ))
    end
    
    # now if collapse flags are true then try to reclassify
    if pargs["biotype"]
      redbiotypes = reducef( reducebiotype, s )           
      s = [s; join(redbiotypes, ':')]
    end
    if pargs["geneid"]
      redgenes = reducef( reducegeneid, s )
      if length(redgenes) >= 2
        if length(unique(redgenes)) == 1
          s[1] = replace(s[1], r"[PI]", "S") #re-classify
        end
      end
      s = [s; join(redgenes, ':')]
    end
    # # # print to STDOUT # # # 
    println( join( s, '\t' ) )
    # # # # # # # # # # # # # # 
  end
end

function masksplit( str, char )
  rep = replace( str, r"[a-z]+", "" )
  split( rep, char )
end

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ### MAIN  ### ### ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #
function main()
  pargs = parse_cmd()

  # simplify ridiculous split string type
  # in case julia changes by version
  typealias SubType typeof( split("a:b", ':') )

  djunc  = Dict{ASCIIString,Bool}()
  dstore = Dict{Char,Array{SubType,1}}() 
  dclass = pargs["load"] == nothing ? Dict{ASCIIString,Char}() : loaddict(pargs["load"])
  dbarcode = pargs["barcode"] == nothing ? Dict{ASCIIString,ASCIIString}() : loaddict(pargs["barcode"])

  const seqInd = 11 # sequence index of .lig
  const geneInd = 4 # gene-id index
  const nameInd = 10 # readName index
  const ligqInd = 13 # LIGQ index

  pargs["load"] == nothing && (sizehint!(dclass, 1000000))
  # first iteration through file, store data, set structures
  for i::ASCIIString in eachline( STDIN )
    s = split(chomp(i), '\t')
    #genes = split(s[geneInd], ':')
    genes = reducef( reducegeneid, s ) # TODO: this gets called twice. fix that.

    # test if this is a unique junction/readset, short circuit by --uniq
    if !pargs["uniq"] || isUniqueJunc!( djunc, s[seqInd], genes, dbarcode, s[nameInd] )
      #set data
      @assert( length(s[1]) == 1 )
      curclass = s[1][1]
      #check ligq for intramolecular matches
      if curclass != 'S' && ismatch(s[ligqInd], r"\[\d+\]")
        s[1] = replace(s[1], r"[PI]", "S")
        curclass = 'S'
      end
      dush!( dstore, curclass, s, arrType=typeof(s) )
      
      # check if we have already loaded "dclass"
      seqs = masksplit(s[seqInd], '_')
      if length(seqs) > 1 && pargs["load"] == nothing
        setClass!( dclass, curclass, seqs, size=pargs["size"], set=true )
      end
    end
  end

  if pargs["save"] != nothing
#    save(pargs["save"], "dclass", dclass)
    savedict( pargs["save"], dclass )
    quit()
  end

  # finish and print.
  reclassReduceAndPrint( dclass, dstore, pargs, seqInd )
end
#########
main()
########
