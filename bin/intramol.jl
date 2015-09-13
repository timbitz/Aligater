#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 6/23/2015
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
function read_lig_parse_reg( io )
   for l::ASCIIString in eachline( io )
      s = split(l, '\t')
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
