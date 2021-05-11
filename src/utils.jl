
"""
Taken from StatisticalTraits.jl
"""
function is_uppercase(char::Char)
    i = Int(char)
    i > 64 && i < 91
end

"""
    snakecase(str, del='_')

Taken from StatisticalTraits.jl
Return the snake case version of the abstract string or symbol, `str`, as in

    snakecase("TheLASERBeam") == "the_laser_beam"

"""
function snakecase(str::AbstractString; delim='_')
    snake = Char[]
    n = length(str)
    for i in eachindex(str)
        char = str[i]
        if is_uppercase(char)
            if i != 1 && i < n &&
                !(is_uppercase(str[i + 1]) && is_uppercase(str[i - 1]))
                push!(snake, delim)
            end
            push!(snake, lowercase(char))
        else
            push!(snake, char)
        end
    end
    return join(snake)
end

snakecase(s::Symbol) = Symbol(snakecase(string(s)))


"""
Taken from MLJBase.jl
"""
function generate_name!(M, existing_names)
    str = split(string(M), '{') |> first
    candidate = split(str, '.') |> last |> snakecase |> Symbol
    candidate in existing_names ||
        (push!(existing_names, candidate); return candidate)
    n = 2
    new_candidate = candidate
    while true
        new_candidate = string(candidate, n) |> Symbol
        new_candidate in existing_names || break
        n += 1
    end
    push!(existing_names, new_candidate)
    return new_candidate
end