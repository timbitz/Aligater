#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 9/13/2015
   e-mail: tim.sterne.weiler@utoronto.ca
=#

using ArgParse

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--bed"
      help = "toggle output as bed format (default on)"
      arg_type = :store_true
  end
  return parse_args(s)
end


###################################################################
function splitgb( coord )
  spl = split( coord, ':' )
  @assert( length(spl) > 1 )
  nums = map( x->parse(Int, x), split( spl[2], '-' ) )
  length(spl) >= 3 ? (spl[1], nums..., spl[3]) : (spl[1], nums...)
end

function read_lig_parse_reg( io )
   const geneind = 25
   const posind = 17
   for l::ASCIIString in eachline( io )
      s = split(l, '\t')
      (g1,g2) = split( s[geneind], ':' )
      g1 == g2 || continue
      (c1,c2) = map( x->splitgb(x), split( s[posind], ',' ) )
      # unfinished as of 10/2015
   end
end

function print_bed( io, reg )
   for bedtup in reg
      bedstr = join(bedtup, '\t')
      println( io , bedstr )
   end
end

function main()
   pargs = parse_cmd()

   reg = read_lig_parse_reg( STDIN )
   if pargs["bed"]
      print_bed( STDOUT, reg )
   end
end

main()
