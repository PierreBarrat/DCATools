using Test
using Crayons.Box


fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
anyerrors = false

my_tests = [
	"IO/test.jl",
	"modeltools/test.jl",
	"aligntools/test.jl",
	"MCMC/test.jl",
]
println("Running tests:")

for my_test in my_tests
    	println("\t", BOLD(BLUE_FG("TESTING")), " $(my_test)")
        include(my_test)
    # try
    # catch e
    #     global anyerrors = true
    #     println("\t\033[1m\033[31mFAILED\033[0m: $(my_test)")
    #     if fatalerrors
    #         rethrow(e)
    #     elseif !quiet
    #         showerror(stdout, e, backtrace())
    #         println()
    #     end
    # end
    println()
end

if anyerrors
    throw("Tests failed")
end
