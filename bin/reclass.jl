#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 3/16/2015
   e-mail: tim.sterne.weiler@utoronto.ca
 
   dev on stable julia 0.3.8
=# 

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ## INITIALIZATION  ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

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
  end
  return parse_args(s)
end

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ###  FUNCTIONS ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

include("dictext.jl") #--> dush!, dinc!, dnorm!

function isUniqueJunc!( used::Dict{ASCIIString,Bool}, seq, genes )
  m = match(r"([AGCTUN]{5}_[AGCTUN]{5})", seq)
  cap = m.captures[1]
  m = match(r"([AGCTUN]+)_([AGCTUN]+)", seq)
  lLen,rLen = length(m.captures[1]), length(m.captures[2])
  geneid = sort( genes )
  key = join(geneid, ":") * "_$cap:$lLen:$rLen"
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
function setClass!( doubleDict::Dict{ASCIIString,Dict{ASCIIString,Char}}, class::Char, seqs; size=32 )
  @assert( length(seqs) >= 2 )
  classHeirPair( a::Char, b::Char ) = a < b ? b : a  #--> Char
  classHeirArray( arr::Array{Char,1} ) = shift!( reverse( sort( arr ) ) ) #--> Char
  retarray = Char[]
  lenA,lenB = map(length, seqs)
  # pad a sequence if it is shorter than our target window size
  seqA = seqs[1] * repeat(".", max( size - lenA, 0 ))
  seqB = seqs[2] * repeat(".", max( size - lenB, 0 ))
  # reset the lengths after padding
  lenA,lenB = length(seqA),length(seqB)
  it=max(1, floor(size/2)) # calculate the step size, as half the window size
  # iterate through the cartesian product of windows in seqA and seqB by it step size
  for i = 1:it:(lenA-size)+1, j = 1:it:(lenB-size)+1
    winA = seqA[i:i+size-1] # access substrings 
    winB = seqB[j:j+size-1]
    keyA,keyB = winA < winB ? (winA,winB) : (winB,winA)
    # use non-exported ht_keyindex function to avoid over indexing hash
    indxA = Base.ht_keyindex( doubleDict, keyA )
    if indxA <= 0 #initialize hash if necessary
      doubleDict[keyA] = Dict{ASCIIString,Char}()
      indxA = Base.ht_keyindex( doubleDict, keyA )
    end
    dictA = doubleDict.vals[indxA]
    indxB = Base.ht_keyindex( dictA, keyB )
    if indxB <= 0
      dictA[keyB] = class #set val if none exists
      length(retarray) <= 0 && push!(retarray, class)
    else
      curVal = classHeirPair( dictA.vals[indxB], class )
      dictA.vals[indxB] = curVal
      push!(retarray, curVal)
    end
  end
  @assert( length(retarray) > 0 )
  length(retarray) > 1 ? classHeirArray( retarray ) : retarray[1]
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
  ismatch(r"(srp|sn|t|r)RNA", biotype) && return(biotype)
  ismatch(r"SNOR|SCAR", geneid) && return("snoRNA")
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
  const geneInd = 3 # gene-id index
  const biotInd = 6 # biotype index
  const repnInd = 7 # repeat name index
  const repcInd = 8 # repeat class index
  
  function safesplit( arr, ind; char=':' )
    @assert( length(arr) >= ind )
    split( arr[ind], char )
  end

  gspl = safesplit( s, geneInd )
  bspl = safesplit( s, biotInd )
  nspl = safesplit( s, repnInd )
  cspl = safesplit( s, repcInd )

  lenLim = min( map(length, (gspl, bspl, nspl, cspl)) )
  res = ASCIIString[]
  for i in 1:lenLim
    push!(res, func(gspl[i], bspl[i], nspl[i], cspl[i]))
  end
  res
end #--> ASCIIString[]

function reclassReduceAndPrint( dclass::Dict, dstore::Dict, pargs, seqInd )
  for class in ['S','P','A','I'], s in dstore[class]
    seqs = split(s[seqInd], '_')
    if length(seqs) > 1 # try to reclassify
      s[1] = setClass!( dclass, 'A', seqs )
    end
    
    # now if collapse flags are true then try to reclassify
    if pargs["biotype"]
      redbiotypes = join(reducef( reducebiotype, s ), ':')
      s = [s, redbiotypes]
    end
    if pargs["geneid"]
      redgenes = join(reducef( reducegene, s ), ':')
      s = [s, redgenes]
    end
    # # # print to STDOUT # # # 
    println( join( s, '\t' ) )
    # # # # # # # # # # # # # # 
  end
end

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ### MAIN  ### ### ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #
function main()
  pargs = parse_cmd()

  # simplify ridiculous split string type
  # in case julia changes by version
  const stype = typeof( split("a:b", ':') )

  djunc = Dict{ASCIIString,Bool}()
  dstore = Dict{Char,Array{stype,1}}()
  dclass = Dict{ASCIIString,Dict{ASCIIString,Char}}()

  const seqInd = 11 # sequence index of .lig
  const geneInd = 3 # gene-id index

  # first iteration through file, store data, set structures
  for i in eachline( STDIN )
    s = split(chomp(i), '\t')
    genes = split(s[geneInd], ':')

    # test if this is a unique junction/readset
    if !pargs["uniq"] || isUniqueJunc!( djunc, s[seqInd], genes )
      #set data
      @assert( length(s[1]) == 1 )
      curclass = s[1][1]
      dush!( dstore, curclass, s, arrType=typeof(s) )
    
      seqs = split(s[seqInd], '_')
      if length(seqs) > 1 
        setClass!( dclass, curclass, seqs )
      end
    end
  end

  # finish and print.
  reclassReduceAndPrint( dclass, dstore, pargs, seqInd )
end
#########
main()
#########
