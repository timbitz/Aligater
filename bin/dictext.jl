#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 3/16/2015
   e-mail: tim.sterne.weiler@utoronto.ca

=#

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

