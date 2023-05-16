# circuitstring conversion to and from impedance.py

function impy_to_ec(circuitstring::String)
    for (f, t) in zip(["(", ")", "p","o"], ["[",  "]", "" ,""])
    circuitstring=replace(circuitstring, f=>t)
    end
    return number_circuit(denumber_circuit(circuitstring))
end

function ec_to_impy(circuitstring::String)
    for (f, t) in zip(["[",  "]", "W" ],["p(", ")", "Wo"])
        circuitstring = replace(circuitstring, f=>t)
        end
        return circuitstring
end


# Importing measurement data from Gamry files.

function ReadGamry(fn_DTA::AbstractString)
    if ~isfile(fn_DTA)
        error("No file found: \"", fn_DTA, "\"!")
    end
    txt_lines = readlines(fn_data, keep=true)
    idx_header_line = nothing
    char_dicimal_delim = '.'
    r_search_pattern_first_line = Regex("^ZCURVE")
    r_comma_in_numbers = Regex("\t\\d+,\\d+")
    idx_line = 0
    for i_line in txt_lines
        idx_line += 1
        if occursin(r_search_pattern_first_line, i_line)
            idx_header_line = idx_line + 3
            @info("first line found :-) #: ", idx_header_line)
        end
        if idx_line == idx_header_line
            println("1st data Line: ", i_line)
        end
        if ~isnothing(idx_header_line) && (idx_line == idx_header_line + 1)
            if occursin(r_comma_in_numbers, i_line)
                char_dicimal_delim = ','
                @info("Comma is decimal delimiter!")
            end
        end
    end
    # --- extract data
    if isnothing(idx_header_line)
        error("Begin of data block not found!")
    else
        # --- original: Pt          Time	Freq	 Zreal	   Zimag	 Zsig	      Zmod	   Zphz	    Idc	    Vdc	    IERange
        s_header = ["", "Point",    "Time", "Frequ", "Z_real", "Z_imag", "Amplitude", "Z_mod", "Z_phz", "I_dc", "V_dc", "IE_range"]
        filecontent = CSV.File(fn_data; 
            header = s_header,
            skipto = idx_header_line, 
            decimal = char_dicimal_delim, 
            delim = "\t");
    end
    # --- check column missmatch:
    n_colums = size(filecontent.names)[1]
    if n_colums != size(s_header)[1]
        error("Column Missmatch")
    end
    _frequ_data = filecontent[:Frequ]
    _Z_data = ComplexF64.(filecontent[:Z_real], filecontent[:Z_imag])
    idx_sorted = sortperm(_frequ_data)
    # ---
    return _frequ_data[idx_sorted], _Z_data[idx_sorted]
end

# Suppressing undesired warning messages

macro suppress(block)
    quote
        if ccall(:jl_generating_output, Cint, ()) == 0
            original_stdout = stdout
            out_rd, out_wr = redirect_stdout()
            out_reader = @async read(out_rd, String)

            original_stderr = stderr
            err_rd, err_wr = redirect_stderr()
            err_reader = @async read(err_rd, String)

            # approach adapted from https://github.com/JuliaLang/IJulia.jl/pull/667/files
            logstate = Base.CoreLogging._global_logstate
            logger = logstate.logger
            if :stream in propertynames(logger) && logger.stream == original_stderr
                new_logstate = Base.CoreLogging.LogState(typeof(logger)(err_wr, logger.min_level))
                Core.eval(Base.CoreLogging, Expr(:(=), :(_global_logstate), new_logstate))
            end
        end

        try
            $(esc(block))
        finally
            if ccall(:jl_generating_output, Cint, ()) == 0
                redirect_stdout(original_stdout)
                close(out_wr)

                redirect_stderr(original_stderr)
                close(err_wr)

                if :stream in propertynames(logger) && logger.stream == stderr
                    Core.eval(Base.CoreLogging, Expr(:(=), :(_global_logstate), logstate))
                end
            end
        end
    end
end