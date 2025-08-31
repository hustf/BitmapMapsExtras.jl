"""
# Example

```
julia> Domain(1.0, 2.0, 3.0, 4.0)
Domain(1.0, 2.0, 3.0, 4.0)

julia> @define_glyphspec_show Domain

julia> Domain(1.0, 2.0, 3.0, 4.0)
Domain(minx=1.0, miny=2.0, maxx=3.0, maxy=4.0)
```
"""
macro define_show_with_fieldnames(T)
    quote
        function Base.show(io::IO, ::MIME"text/plain", x::$(esc(T)))
            print(io, $(string(T)) * "(")
            fields = fieldnames($(esc(T)))
            for (i, field) in enumerate(fields)
                print(io, string(field) * "=")
                show(io, getfield(x, field))
                if i < length(fields)
                    print(io, ", ")
                end
            end
            print(io, ")")
        end
    end
end