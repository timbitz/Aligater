#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 3/16/2015
   e-mail: tim.sterne.weiler@utoronto.ca

   Calculate enrichment of chimeras over expression level background
   using multinomial draws of n=2 to calculate a second multinomial
   distribution of gene pairs which is used to compare observed frequency
   to expected where each pair is binomially distributed 
   k occurences of chimera (i,j) is k[(i,j)]~ Bin( n, p[(i,j)] ) 
=# 

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ## INITIALIZATION  ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

println(STDERR, "[aligater enrich] Loading Packages..")

using ArgParse
using StatsBase
using Distributions

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--back", "-b"
      help = "lig file to build background from"
      arg_type = ASCIIString
    "--fore", "-f"
      help = "lig file to evaluate"
      arg_type = ASCIIString
    "--gi"
      help = "index for geneid:geneid field in foreground"
      arg_type = Int
      default = 25
    "--ci"
      help = "index for counting field if applicable (non integer fields increment by 1 instead)"
      arg_type = Int
    "--bi"
      help = "index for background calibration (default is [end] of collapsed lig without any post)"
      arg_type = Int
      default = 19
    "--nd"
      help = "comma delmited set of strings to substitute for wildcard.  The ratio of each pair is output"
      arg_type = ASCIIString
    "--filt"
      help = "pattern filter [syntax= column:regex]"
      arg_type = ASCIIString
      default = "1:."
    "--alpha"
      help = "maximum p-value for output, only applied to the numerator (first element) of --nd if supplied"
      arg_type = Float64
      default = 1.0
    "--nc"
      help = "normalization constants, applied to each of --nd as a comma separated list."
      arg_type = ASCIIString
  end
  return parse_args(s)
end

macro vprint(msg)
# todo
end

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ###  FUNCTIONS ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

function parse_colfilt(toparse::ASCIIString)
  cap = match(r"(\d+)\:(.*)", toparse).captures
  @assert( length(cap) == 2 )
  col = int(cap[1])
  reg = Regex(cap[2])
  col, reg
end #--> (Integer, Regex)

# dictionary value increment
# this should be standard imo.
function dinc!(dict::Dict, key, val=1) 
  if !haskey(dict, key)
    dict[key] = val
  else
    @assert( isa(dict[key], Number) )
    dict[key] += val
  end
end #--> nothing

# dictionary [] create or push if exists
# this also should be standard imo
function dush!(dict::Dict, key, val; arrType = Any)
  @assert( isa(val, arrType) )
  if !haskey(dict, key)
    dict[key] = arrType[ val ]
  else
    push!(dict[key], val)
  end
end #--> nothing


# normalize numeric dictionary values
function dnorm!(dict::Dict)
  const n = sum( collect( values(dict) ) )
  @assert( n > 0 )
  for i in keys(dict)
    @assert( isa(dict[i], Number) )
    dict[i] /= n 
  end #normalize
end #--> nothing

#= make joint interaction distribution 
     P(x,y) = [ x * y  when x != y
              |
              [   0    when x == y 
=#
## Use this function to find the denominator for x != y
#  By dividing each x * y pair by 1-Sum(x == y)
function sumSamePairs(pset::Dict)
  sum = 1.0
  for i in keys(pset)
    sum -= 2(pset[i] ^ 2)
  end
  sum
end #--> Float64

# make a pseudoCount
function pseudoCnt(pset::Dict)
  pseudo = minimum( collect( values( pset ) ) )^2
  pseudo * 0.5
end #--> Float

# get probability of pairs
function jointInteracProb(pset::Dict, key1, key2, denom = 1.0, pseudo = 0.0)
  if(!haskey(pset, key1) || !haskey(pset, key2))
    return( pseudo )
  end
  2(pset[key1] * pset[key2]) / denom
end #--> Float

function calculateBinomialStats(forecnt::Dict, backprob::Dict)
  const n = sum( collect( values(forecnt) ) )
  const den = sumSamePairs(backprob)
  const pseudo = pseudoCnt(backprob)
  const bonfCor::Int = length(keys(forecnt))

  foreprob = Dict{Tuple,Float64}()
  retval = Dict{(String,String),Tuple}()
 
  obsDen = 0.0

  for (k1,k2) in keys(forecnt)
    prob = jointInteracProb(backprob, k1, k2, den, pseudo)
    obsDen += prob
    foreprob[(k1,k2)] = prob
  end

  for (k1,k2) in keys(forecnt)
    prob = foreprob[(k1,k2)] / obsDen
    bin = Binomial(int(n), prob)
    k = forecnt[(k1,k2)]
    pval = ccdf(bin, k-1)
    pval *= bonfCor
    pval = (pval > 1.0) ? one(pval) : pval
    expval = prob * n
    retval[(k1,k2)] = (k, expval, pval, k/expval)
  end
  retval
end #--> Dict{(String,String), Tuple}

# if the pairs 
function loadInteractionFile(io::IOStream, gInd::Int, cntInd::Int, col, reg)
  cset = Dict{Tuple,Float64}()
  sizehint!(cset, 100000) # allocate 100K pairs
  for l::ASCIIString in eachline( io )
    s = split(chomp(l), '\t')
    if isa(col, Integer) && isa(reg, Regex) && length(s) >= col
      ismatch(reg, s[col]) || continue
    end
    # make sure we don't get out of bounds error
    (ismatch(r":", s[gInd])) || continue
    k1,k2 = sort( split(s[gInd], ':') )
    cnum = parse( s[cntInd], raise=false )
    dinc!(cset, (k1,k2), isa(cnum, Number) ? convert(Float64, cnum) : 1.0)
  end
  cset
end #--> Dict{Tuple,Float64}

# this does io and creates a cset structure of type Dict
function loadBackgroundFile(io::IOStream, ind::Int; keyType=String, valType=Float64)
  cset = Dict{keyType,valType}()
  sizehint!(cset, 20000) # approx number of genes
  for i::ASCIIString in eachline(io)
    s = split(chomp(i), '\t')
    #=genes = replace(s[ind], ":*NA:*", "") |> x->split(x, ':')
    genes = split(s[ind], ':')
    length(genes) > 0 && dinc!(cset, genes[1], one(valType)) =# # 04/15/15 deprecate for performance
    dinc!(cset, s[ind], one(valType))
  end
  cset
end #--> Dict{keyType,valType}

# This takes a cmd line string that contains pairs=> "column:type,column:type"
# where column is an int and type is a Char of [sfdn] 
function varStatSummary( io, refInd::Int, varstr::ASCIIString, col, reg )
  # helper functions:
  function parseVarStr( str::ASCIIString )
    sets = split(str, ',')
    ret  = (Int,Char)[]
    for i in sets
      col,typ = map(parse, split(i, ':'))
      @assert( ismatch( r"[cf]", typ[1] ) )
      push!(ret, (col, typ[1]) )
    end
    ret
  end #--> Array{(Int,Char),1}

  function letterToType( letter::Char )
    letter == 'p' && return (ASCIIString, ASCIIString)
    letter == 's' && return ASCIIString
    letter == 'f' && return Float64
    letter == 'd' && return Int64
    letter == 'n' && return Number
    Any
  end #--> Type{T}
  
  isOrdered( t::Tuple ) = length(t) <= 1 || t[1] <= t[2] ? true : false

  # varStatSummary code:
  varSet = parseVarStr( varstr )
  summary = Dict{Int,Dict}()
  for l::ASCIIString in eachline(io)
    s = split(chomp(l), '\t')
    # follow the same filter rules as interactionFile:
    if isa(col, Integer) && isa(reg, Regex) && length(s) >= col
      ismatch(reg, s[col]) || continue
    end
    # get reference key.
    refKey = sort( split(s[refInd], ':') )
    k = tuple(refKey...) # main refKey k
    # iterate through column & type pairs
    for (i,c) in varSet
      cType = letterToType( c ) # convert char to type
      parSi = parse(s[i])  # parse value initially
      cVal  = cType <: Tuple ? begin # if tuple check if properly ordered
                                 (a,b) = split(parSi, ':')
                                  isOrdered( k ) ? (a,b) : (b,a)
                               end : parSi #otherwise initial parse was fine
      @assert( isa(cVal, cType) )
      if !haskey(summary, i) 
        vType = cType <: ASCIIString || cType <: Tuple ? cType : Array{cType,1}
        summary[i] = Dict{typeof(k), vType}()
      end
      if cType <: ASCIIString # set to single value unless already set
        !haskey( summary[i], k ) && (summary[i][k] = cVal)
      else # numeric, so push the current element
        push!( summary[i][k], cVal )
      end
    end 
  end
  summary
end #--> Dict{Int,Dict{Tuple, ? }}

# central function for processing a foreground background filename pair
function loadFilesAndCalculate(forefile::ASCIIString, backfile::ASCIIString, pargs)
  # fetch options.
  const backInd = pargs["bi"]
  const geneInd = pargs["gi"]
  const cntInd  = pargs["ci"] == nothing ? geneInd : pargs["ci"]

  c,r = parse_colfilt(pargs["filt"])

  print(STDERR, "Loading Background $backfile..\n")
  backhndl = open(backfile, "r")
  cset = loadBackgroundFile(backhndl, backInd, valType=Float64)
  close(backhndl)

  print(STDERR, "Processing probability distribution..\n")
  # process the probability distribution
  pset = deepcopy(cset) # multinomial probabilities
  dnorm!( pset ) # normalize dict

  # iterate through foreground
  print(STDERR, "Loading Foreground $forefile..\n")
  forehndl = open(forefile, "r")
  forecnt  = loadInteractionFile(forehndl, geneInd, cntInd, c, r) 
  close(forehndl)

  # calculate p-values and return hash
  calculateBinomialStats( forecnt, pset ) # return stat hash
end
#--> Dict{(String,String), Tuple}

# final print function for stat array
function printStats( io, statarr::Array; normarr=[1,1], alpha=1.0, vardict=Dict() )
  aHsh = statarr[1]
  @assert( isa(aHsh, Dict) )

  # internal print function for variable keys
  function printVardict( io, dict, refkey )
    for col in keys(dict)
      if haskey( dict[col], refkey )
        val = dict[col][refkey]
        if isa( val, Array )
          med = median( val )
          mad = mad( val )
          @printf( io, "\t%d\t%.2f±%.2f", col, med, mad)
        elseif isa( val, Tuple )
          @printf( io, "\t%d\t%s", col, join( val, ',' ) )
        else
          @printf( io, "\t%d\t%s", col, val )
        end
      else
        @printf( io, "\t%d\tNA", col )
      end   
    end
    @printf( io, "\n" )
  end #--> nothing

  if length(statarr) == 1
    # just a single set print output normally
    for (k1,k2) in keys( aHsh )
      k, _, pval, obsExp = aHsh[(k1,k2)]
      aPval <= alpha || continue  #filter p-val
      @printf( io, "%s,%s\t%d\t%.2e\t%.2e", k1, k2, k, pval, obsExp )
      printVardict( io, vardict, (k1,k2) )
    end
  elseif length(statarr) == 2
    bHsh = statarr[2]
    @assert( isa(bHsh, Dict) )
    ksone  = collect( keys( aHsh ) )
    kstwo  = collect( keys( bHsh ) )
    intset = intersect( ksone, kstwo )
    # iterat ehrough shared set.
    for (k1,k2) in ksone
      aK, aExp, aPval, aObsExp = aHsh[(k1,k2)]
      bK, bExp, bPval, bObsExp = haskey(bHsh,(k1,k2)) ? bHsh[(k1,k2)] : (0.5, 0.5, 1.0, 1.0)
      a_b = (aK / normarr[1]) / (bK / normarr[2])
      exp_a_b = (aExp / normarr[1]) / (bExp / normarr[2])
      full_a_b = a_b / exp_a_b
      aPval <= alpha || continue # filter if p-value is above cutoff
      @printf( io, "%s,%s\t%.1f\t%.1f\t%.1f\t%d\t%d\t%.2e\t%.2e", k1, k2, full_a_b, a_b, exp_a_b, aK, bK, aPval, bPval )
      printVardict( io, vardict, (k1,k2) )
    end #endfor
  end #endif
end #--> nothing

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ### MAIN  ### ### ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #
function main()
  pargs = parse_cmd()

  # check options.
  ndArray = pargs["nd"] == nothing ? [""] : split(pargs["nd"], ",")
  normArray = pargs["nc"] == nothing ? [1,1] : map(parse, split(pargs["nc"], ","))

  stats = Dict[]
  
  varstat = nothing  #scope

  # iterate through delimited names list
  for (ndIter,nd) in enumerate(ndArray)

    # if wilcard present, substitute with current pattern
    backfile = replace(pargs["back"], "%", nd)
    forefile = replace(pargs["fore"], "%", nd)

    if ndIter == 1
      fh = open(forefile, 'r')
      varstat = varStatSummary( fh, pargs["gi"], parse_colfilt(pargs["filt"]) )
      close( fh )
    end

    # calculate p-values
    stathsh = loadFilesAndCalculate( forefile, backfile, pargs )
    
    push!(stats, stathsh) 
  end  

  @assert(0 < length(stats) <= 2, "--nd does not contain properly formatted arguments!")
  printStats( STDOUT, stats, alpha=pargs["alpha"], normarr=normArray, vardict=varstat )

end
#####################################################
# main execution here
main()
# eof