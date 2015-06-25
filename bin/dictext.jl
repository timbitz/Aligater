#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 3/16/2015
   e-mail: tim.sterne.weiler@utoronto.ca

=#
using GZip

# slow performance of HDF5 has promted [savedict, loaddict]:
function savedict{K,V}( file::String, dict::Dict{K,V} )
  GZip.open(file, "w") do fh
    println(fh, "@Dict{" * string(K) * "," * string(V) * "}:" * string(length(dict)))
    for k in keys(dict)
      println(fh, string(k) * "˧" * string(dict[k]))
    end
  end
end #--> nothing

function loaddict(file::String)
   myconvert(x::Type{Char}, y) = string(y)[1]
   myconvert(x::Type{ASCIIString}, y) = convert(ASCIIString, string(y))
   myconvert{T <: Number, S <: String}(x::Type{T}, y::S) = parse(T, y)
   myconvert(x, y) = convert(x, y)
   #hint = parse(chomp(readall(pipe(`zcat $file`, `wc -l`)))) # --deprecated 6/10/15 tsw
   dict = nothing
   GZip.open(file, "r") do fh
    heads = split(chomp( readline(fh) ), r"@Dict{|,|}:")
    @assert( length(heads) == 4 )
    @assert( ismatch(r"^[A-Z|a-z|0-9]+$", heads[2]) ) # must look like a type
    @assert( ismatch(r"^[A-Z|a-z|0-9]+$", heads[3]) ) # before we allow parse()
    ktype = eval(parse(heads[2]))
    vtype = eval(parse(heads[3]))
    dict = Dict{ktype,vtype}()
    hint  = parse(Int, heads[4])
    sizehint!(dict, hint)
    for l in eachline(fh)
      sub = split(chomp(l), '˧')
      key = myconvert(ktype, sub[1])
      val = myconvert(vtype, sub[2])
      dict[key] = val
    end
  end
  dict
end #--> Dict{K,V}

# dictionary value increment
# this should be standard imo.
function dinc!{K,V}(dict::Dict{K,V}, key::K, val::V=one(V))
  indx = Base.ht_keyindex(dict, key)
  if indx <= 0
    dict[key] = val
  else
    @assert( isa(dict.vals[indx], Number) )
    dict.vals[indx] += val
  end
end #--> nothing

# dictionary [] create or push if exists
# this also should be standard imo
function dush!{K,V}(dict::Dict{K,V}, key::K, val; arrType = Any)
  indx = Base.ht_keyindex(dict, key)
  @assert( isa(val, arrType) )
  if indx <= 0
    dict[key] = arrType[ val ]
  else
    push!(dict.vals[indx], val)
  end
end #--> nothing


# normalize numeric dictionary values
function dnorm!{K, V<:Number}(dict::Dict{K,V})
  const n = sum( collect( values(dict) ) )
  @assert( n > 0 )
  for i in keys(dict)
    dict[i] /= n
  end #normalize
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

