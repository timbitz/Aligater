#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 3/16/2015
   e-mail: tim.sterne.weiler@utoronto.ca

=# 

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ## INITIALIZATION  ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

global head = "[aligater table]:"
println(STDERR, "$head Loading Packages..")

using ArgParse
using StatsBase

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--gi"
      help = "index for geneid:geneid field in foreground"
      arg_type = Int
      default = 1
    "--ci"
      help = "index for counting field if applicable (non integer fields increment by 1 instead)"
      arg_type = Int
      default = 2
    "--filt","-f"
      help = "pattern filter [syntax= column:regex]"
      arg_type = ASCIIString
      default = "1:."
    "--alpha","-a"
      help = "maximum p-value for output, only applied to the numerator (first element) of --nd if supplied"
      arg_type = Float64
      default = 1.0
    "--delim","-d"
      help = "delimiter for paired index entry"
      arg_type = ASCIIString
      default = ":"
    "--null","-n"
      help = "null value to insert for missing pairs"
      arg_type = ASCIIString
      default = "X"
    "--outsep","-s"
      help = "output separator for table"
      arg_type = ASCIIString
      default = ","
    "--integer","-i"
      help = "only count integers, try to convert if --ci is not an integer"
      action = :store_true
    "--norepeat","-r"
      help = "split names and remove _\\S+"
      action = :store_true
    "--sumOnly"
      action = :store_true
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
end #--> Tuple{Integer, Regex}

# if the pairs --TODODOC
function load_interactionfile(io, gInd::Int, cntInd::Int, delim, col, reg; ctype=Float64, spltflag=false)
  cset = Dict{AbstractString,Dict{AbstractString,ctype}}()
  #sizehint!(cset, 100)
  sum = zero(ctype)
  for l in eachline( io )
    s = split(chomp(l), '\t')
    if isa(col, Integer) && isa(reg, Regex) && length(s) >= col
      ismatch(reg, s[col]) || continue
    end
    # make sure we don't get out of bounds error
    (ismatch(Regex(delim), s[gInd])) || continue
    k1,k2 = spltflag ? map( x->split(x, '_')[1], sort( split(s[gInd], delim) ) ) : sort( split(s[gInd], delim) )

    if ismatch(r"sense|transcript|protein|mRNA|sens|int|lncRNA", k1) || ismatch(r"snoRNA|miRNA", k2)
      k1,k2 = k2,k1
    end
 
    cnum = parse( s[cntInd], raise=false )
    if !haskey(cset, k1)
      cset[k1] = Dict{AbstractString,ctype}()
    end
    indict = cset[k1]
    toinc = isa(cnum, Number) ? convert(ctype, cnum) : one(ctype)
    sum += toinc
    dinc!(indict, k2, toinc)
  end
  cset,sum
end #--> Dict{ASCIIString,Dict{ASCIIString,ctype}},sum::ctype

# tested.
function print_table{K,V}( io, table::Dict{K,Dict{K,V}}; delim=",", null="X" )
 allk = Set{K}()
 for i::K in keys(table), j::K in keys(table[i]) 
   push!(allk, i)
   push!(allk, j)
 end
 header = "TAB" * delim * join( allk, delim )
 println( header )
 for i::K in allk
   print( io, i )
   if haskey( table, i )
     indict = table[i]
     for j::K in allk
       val = haskey(indict, j) ? string(indict[j]) : null
       print( io, delim * val ) 
     end
   else 
     print( io, repeat( delim * null, length(allk) ) )
   end #if
   print( io, "\n" )
 end
end

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ### MAIN  ### ### ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #
function main()
  pargs = parse_cmd()

  col,reg = parse_colfilt( pargs["filt"] ) 
  const mytype = pargs["integer"] ? Int64 : Float64;

  table,sum = load_interactionfile( STDIN, pargs["gi"], pargs["ci"], pargs["delim"], col, reg, ctype=mytype, spltflag=pargs["norepeat"] ) 
  pargs["sumOnly"] && println(STDOUT, "$(sum)")
  
  print_table( STDOUT, table; delim=pargs["outsep"], null=pargs["null"] ) 
end
#####################################################
# main execution here
main()
# eof
