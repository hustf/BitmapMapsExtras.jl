# For testing graphically
# 
using SHA # Test dependency, temporarily added as a direct dependency.

using BitmapMapsExtras: display_if_vscode, Colorant
using BitmapMapsExtras: channelview

function hash_image(img)
    io = IOBuffer()
    # Record minimal metadata
    write(io, UInt64(ndims(img)))
    foreach(s -> write(io, UInt64(s)), size(img))
    tstr = String(nameof(eltype(img)))
    write(io, UInt64(sizeof(codeunits(tstr))))
    write(io, codeunits(tstr))
    # Record pixel bytes in a deterministic order
    a = channelview(img)                    # (channels, H, W) or (channels, …)
    write(io, reinterpret(UInt8, vec(a)))   # linearized, column-major
    bytes2hex(sha1(take!(io)))
end

const COUNT = Ref(0)
(::typeof(COUNT))() = COUNT[] += 1

function is_hash_stored(img, vhash)
    if eltype(img) <: Colorant
        display_if_vscode(img)
    end
    if isempty(vhash) || (length(vhash) < COUNT[])
        push!(vhash, hash_image(img))
        # This is for pasting into the test criterion, for later
        # (provided that the output IS ok!)
        s = "vhash = " * string(vhash)
        printstyled("\n " * s * "\n", color = :176)
        if isinteractive()
            clipboard(s)
        end
        # This ensures the rest of the tests in this set branch here, too.
        # Updating the test results is quick, with one paste operation per testset.
        COUNT(); COUNT()
        return false
    else
        i = COUNT() # Increases COUNT[]
        hsh = hash_image(img) 
        result = (hsh == vhash[i])
        if ! result
            printstyled("\"" * string(hsh) * "\", ", color = :light_red)
        end
        return result
    end
end
nothing
