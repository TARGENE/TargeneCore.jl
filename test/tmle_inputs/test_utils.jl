function cleanup(;prefix="final.")
    for file in readdir()
        if startswith(file, prefix)
            rm(file)
        end
    end
end