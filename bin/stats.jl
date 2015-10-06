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

global head = "[aligater stats]:"
println(STDERR, "$head Loading Packages..")

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
    "--vs"
      help = "extra variable string etc.., 'col:type,col:type', where type = [pcfdn]"
      arg_type = ASCIIString
      default = ""
    "--alpha"
      help = "maximum p-value for output, only applied to the numerator (first element) of --nd if supplied"
      arg_type = Float64
      default = 1.0
    "--nc"
      help = "normalization constants, applied to each of --nd as a comma separated list."
      arg_type = ASCIIString
    "--rpm"
      help = "function to combine the background expression values for RPM column"
      arg_type = ASCIIString
      default = "min"
  end
  return parse_args(s)
end


## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ###  FUNCTIONS ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

include("dictext.jl")

# this accepts a --filt string argument in the form 'column:regexstring'
function parse_colfilt(toparse::ASCIIString)
  cap = match(r"(\d+)\:(.*)", toparse).captures
  @assert( length(cap) == 2 )
  col = parse(Int, cap[1])
  reg = Regex(cap[2])
  col, reg
end #--> (Integer, Regex)

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
function jointInteracProb{K,V}(pset::Dict{K,V}, key1::K, key2::K, denom = 1.0, pseudo = 0.0)
  if(!haskey(pset, key1) || !haskey(pset, key2))
    return( pseudo )
  end
  2(pset[key1] * pset[key2]) / denom
end #--> Float

#= This function takes a foreground count set for gene pairs
   and a background probability set and derives a probability
   set for those pairs and computes binomial p-values for the
   foreground counts =#
function calculateBinomialStats{T <: Tuple, S <: AbstractString}(forecnt::Dict{T,Float64}, backprob::Dict{S,Float64})
  const n = sum( collect( values(forecnt) ) )
  const den = sumSamePairs(backprob)
  const pseudo = pseudoCnt(backprob)
  const bonfCor = length(keys(forecnt))

  foreprob = Dict{Tuple,Float64}()
  retval = Dict{Tuple{ASCIIString,ASCIIString},Any}()
 
  obsDen = 0.0

  for (k1,k2) in keys(forecnt)
    prob = jointInteracProb(backprob, k1, k2, den, pseudo)
    obsDen += prob
    foreprob[(k1,k2)] = prob
  end

  for (k1,k2) in keys(forecnt)
    prob = foreprob[(k1,k2)] / obsDen
    bin = Binomial(Int(n), prob)
    k = forecnt[(k1,k2)]
    pval = ccdf(bin, k-1)
    pval *= bonfCor
    pval = (pval > 1.0) ? one(pval) : pval
    expval = prob * n
    retval[(k1,k2)] = (k, expval, pval, k/expval)
  end
  retval
end #--> Dict{(String,String), Any}

# if the pairs --TODODOC
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
function loadBackgroundFile(io::IOStream, ind::Int; keyType=AbstractString, valType=Float64)
  cset = Dict{keyType,valType}()
  sizehint!(cset, 30000) # approx number of genes
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
function varStatSummary( openwhat, refInd::Int, varstr::ASCIIString, col, reg )
  # helper functions:
  function parseVarStr( str::ASCIIString )
    sets = split(str, ',')
    ret  = Tuple{Int,Char}[]
    for i in sets
      vcol,typ = map(parse, split(i, ':'))
      tchar = string(typ)[1]
      @assert( ismatch( r"[psfdn]", string(tchar) ) )
      push!(ret, (vcol, tchar) )
    end
    ret
  end #--> Array{(Int,Char),1}

  function letterToType( letter::Char )
    letter == 'p' && return typeof(("",""))
    letter == 's' && return ASCIIString
    letter == 'f' && return Float64
    letter == 'd' && return Int64
    letter == 'n' && return Number
    Any
  end #--> Type{T}
  
  reOrder( t::Tuple ) = (length(t) <= 1 || t[1] <= t[2]) ? t : (t[2],t[1])
  toStrings( arr ) = map( x->convert(ASCIIString,x), arr )

  # varStatSummary code:
  varSet = parseVarStr( varstr )
  summary = Dict{Int,Dict}()
  open( openwhat , "r" ) do fh
    for l::ASCIIString in eachline(fh)
      s = split(chomp(l), "\t")
      # follow the same filter rules as interactionFile:
      if isa(col, Integer) && isa(reg, Regex) && length(s) >= col
        ismatch(reg, s[col]) || continue
      end
      # get reference key.
      refKey = sort( split(s[refInd], ':') )
      k = tuple(refKey...)[1:2] # main refKey k
      # iterate through column & type pairs
      for (i,c) in varSet
        cType = letterToType( c ) # convert char to type
        parSi = cType <: Number ? parse(cType, string(s[i]) ) : string( s[i] )
        cVal  = cType <: Tuple ? begin # if tuple check if properly ordered
                                    a,b = split(parSi, ':') |> toStrings
                                    reOrder( (a,b) )
                                 end : convert(cType, parSi)
                                  #otherwise initial parse was fine
        @assert( isa(cVal, cType) )
        if !haskey(summary, i) 
          vType = cType <: ASCIIString || cType <: Tuple ? cType : Array{cType,1}
          summary[i] = Dict{typeof(k), vType}()
        end
        if cType <: ASCIIString || cType <: Tuple # set to single value unless already set
          !haskey( summary[i], k ) && ( summary[i][k] = cVal )
        else # numeric, so push the current element
          dush!( summary[i], k, cVal, arrType=cType )
        end
      end 
    end
  end
  summary
end #--> Dict{Int,Dict{Tuple, ? }}

# central function for processing a foreground background filename pair
#=function loadFilesAndCalculate(forefile::ASCIIString, backfile::ASCIIString, pargs)
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

  # if --rpm then print abundance levels also
  @assert( ismatch(r"^[A-Z|a-z|0-9]+$", pargs["rpm"]) )
  rpmFunc = eval(parse(pargs["rpm"]))
  @assert( typeof(rpmFunc) <: Function )
  rset = deepcopy(pset)
  dnorm!(rset, 1/1000000) # reads per M
  
  # iterate through foreground
  print(STDERR, "Loading Foreground $forefile..\n")
  forehndl = open(forefile, "r")
  forecnt  = loadInteractionFile(forehndl, geneInd, cntInd, c, r) 
  close(forehndl)

  print(STDERR, "Calculating binomial stats..\n")
  # calculate p-values and return hash
  calculateBinomialStats( forecnt, pset, rset, rpmFunc ) # return stat hash
end =#  ##Deprecated 7/16/2015
#--> Dict{(String,String), Any}

# final print function for stat array
function printStats( io, statarr::Array, rpmarr::Array; normarr=[1,1], alpha=1.0, vardict=Dict(), rpmFunc=min )
  aHsh = statarr[1]
  aRpm = rpmarr[1]
  @assert( isa(aHsh, Dict) )
  @assert( isa(aRpm, Dict) )

  # internal print function for variable keys
  function printVardict( io, dict, refkey )
    for col in keys(dict)
      if haskey( dict[col], refkey )
        val = dict[col][refkey]
        if isa( val, Array )
          asfloat = typeof(val) <: Array{Int64,1} ? convert(Array{Float64,1}, val) : val
          aMed = median( asfloat )
          aMad = mad( asfloat )
          @printf( io, "\t%d\t%.2f±%.2f", col, aMed, aMad)
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
      k, kExp, pval, obsExp = aHsh[(k1,k2)]
      aPval <= alpha || continue  #filter p-val
      @printf( io, "%s,%s\t%d\t%.1f\t%.2e\t%.2e", k1, k2, k, k/kExp, pval, obsExp )
      printVardict( io, vardict, (k1,k2) )
    end
  elseif length(statarr) == 2
    bHsh = statarr[2]
    bRpm = rpmarr[2]
    @assert( isa(bHsh, Dict) )
    @assert( isa(bRpm, Dict) )
    ksone  = collect( keys( aHsh ) )
    kstwo  = collect( keys( bHsh ) )
    both = union( ksone, kstwo )
    # function to check existance of and reduce rpm values.
    function apply_rpmfunc{K <: AbstractString,V <: Float64}( hsh1::Dict{K,V}, hsh2::Dict{K,V}, key1::K, key2::K, func::Function)
      begin
        val1 = get(hsh1, key1, zero(V))
        val2 = get(hsh2, key2, zero(V))
        func(val1, val2)
      end
    end #--> V
    # iterat ehrough shared set.
    for (k1,k2) in both
      aK, aExp, aPval, aObsExp = haskey(aHsh,(k1,k2)) ? aHsh[(k1,k2)] : (0.5, 0.5, 1.0, 1.0)
      bK, bExp, bPval, bObsExp = haskey(bHsh,(k1,k2)) ? bHsh[(k1,k2)] : (0.5, 0.5, 1.0, 1.0)
      aRpmVal = apply_rpmfunc( aRpm, aRpm, k1, k2, rpmFunc )
      bRpmVal = apply_rpmfunc( bRpm, bRpm, k1, k2, rpmFunc )
      a_b = (aK / normarr[1]) / (bK / normarr[2])
      exp_a_b = (aExp / normarr[1]) / (bExp / normarr[2])
      full_a_b = a_b / exp_a_b
      aPval <= alpha || continue # filter if p-value is above cutoff
      @printf( io, "%s,%s\t%.1f\t%.1f\t%.1f\t%d\t%.1f\t%d\t%.1f\t%.2e\t%.2e\t%.1f\t%.1f", k1, k2, full_a_b, a_b, exp_a_b, aK, aK/aExp, bK, bK/bExp, aPval, bPval, aRpmVal, bRpmVal )
      printVardict( io, vardict, (k1,k2) )
    end #endfor
  end #endif
end #--> nothing

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ### MAIN  ### ### ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #
function main()
  pargs = parse_cmd()

  const backInd = pargs["bi"]
  const geneInd = pargs["gi"]
  const cntInd  = pargs["ci"] == nothing ? geneInd : pargs["ci"]

  c,r = parse_colfilt(pargs["filt"])

  # check options.
  ndArray = pargs["nd"] == nothing ? [""] : split(pargs["nd"], ",")
  normArray = pargs["nc"] == nothing ? [1,1] : map(parse, split(pargs["nc"], ","))

  @assert( pargs["fore"] != nothing, "$head --fore and --back must be given!")
  @assert( pargs["back"] != nothing, "$head --fore and --back must be given!")

  ndParsed = ASCIIString[]
  stats = Dict[]
  rpms  = Dict[]  

  varstat = nothing  #scope

  # parse --rpm flag to Function
  @assert( ismatch(r"^[A-Z|a-z|0-9]+$", pargs["rpm"]) )
  rpmFunc = eval(parse(pargs["rpm"]))
  @assert( typeof(rpmFunc) <: Function, "typeof(--rpm) must be Julia Function!" )

  # iterate through delimited names list
  for (ndIter,nd) in enumerate(ndArray)

    # if wilcard present, substitute with current pattern
    backfile = replace(pargs["back"], "%", nd)
    forefile = replace(pargs["fore"], "%", nd)

    push!(ndParsed, forefile)

    # calculate p-values
    #stathsh = loadFilesAndCalculate( forefile, backfile, pargs ) ##Deprecated 7/16/2015

    print(STDERR, "Loading Background $backfile..\n")
    backhndl = open(backfile, "r")
    cset = loadBackgroundFile(backhndl, backInd, valType=Float64)
    close(backhndl)

    print(STDERR, "Processing probability distribution..\n")
    # process the probability distribution
    pset = deepcopy(cset) # multinomial probabilities
    dnorm!( pset ) # normalize dict

    # if --rpm then print abundance levels also
    rset = deepcopy(pset)
    dnorm!(rset, 1/1000000) # reads per M

    # iterate through foreground
    print(STDERR, "Loading Foreground $forefile..\n")
    forehndl = open(forefile, "r")
    forecnt  = loadInteractionFile(forehndl, geneInd, cntInd, c, r)
    close(forehndl)

    print(STDERR, "Calculating binomial stats..\n")
    # calculate p-values and return hash
    stathsh = calculateBinomialStats( forecnt, pset ) # return stat hash
    
    push!(stats, stathsh) 
    push!(rpms, rset)
  end  

  if length(pargs["vs"]) > 0
    println(STDERR, "Loading --vs variables")
    c,r = parse_colfilt(pargs["filt"])
    varstat = varStatSummary( `cat $ndParsed` , pargs["gi"], pargs["vs"], c, r )
  end

  @assert(0 < length(stats) <= 2, "--nd does not contain properly formatted arguments!")
  printStats( STDOUT, stats, rpms, alpha=pargs["alpha"], normarr=normArray, vardict=varstat, rpmFunc=rpmFunc )

end
#####################################################
# main execution here
main()
# eof
