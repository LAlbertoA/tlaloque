import REPL
using REPL.TerminalMenus

include("tests_utilities.jl")

tests_path = @__DIR__*"/tests"
if !isdir(tests_path)
    mkdir(tests_path)
end

println("Running test suite for: ")
println("\n\n")

println(" ███████████ █████         █████████   █████          ███████       ██████    █████  █████ ██████████")
println("░█░░░███░░░█░░███         ███░░░░░███ ░░███         ███░░░░░███   ███░░░░███ ░░███  ░░███ ░░███░░░░░█")
println("░   ░███  ░  ░███        ░███    ░███  ░███        ███     ░░███ ███    ░░███ ░███   ░███  ░███  █ ░ ")
println("    ░███     ░███        ░███████████  ░███       ░███      ░███░███     ░███ ░███   ░███  ░██████   ")
println("    ░███     ░███        ░███░░░░░███  ░███       ░███      ░███░███   ██░███ ░███   ░███  ░███░░█   ")
println("    ░███     ░███      █ ░███    ░███  ░███      █░░███     ███ ░░███ ░░████  ░███   ░███  ░███ ░   █")
println("    █████    ███████████ █████   █████ ███████████ ░░░███████░   ░░░██████░██ ░░████████   ██████████")
println("   ░░░░░    ░░░░░░░░░░░ ░░░░░   ░░░░░ ░░░░░░░░░░░    ░░░░░░░       ░░░░░░ ░░   ░░░░░░░░   ░░░░░░░░░░ ")

tests, tests_selected = _run_menu()
println("\n ################################################################################### \n")
println("Running this tests will make a new folder in the test_suit directory called '(test name)_test'.
Inside this folder there will be a Makefile, the user.f90 and parameters.f90 files of the test, and a DATA
folder with the outputs of the test. Also, the test will be compared with expected results and this will
be saved to a file called '(test name)_results.txt'.")
println("\n ################################################################################### \n")

tests_modules = strip.(getindex.(split.(tests[tests_selected], ":"), 1))
tests_names = strip.(getindex.(split.(tests[tests_selected], ":"), 2))
paths_to_tests = joinpath.(tests_path, replace.(tests_names, " " => "_"))

for i in 1:length(tests_names)
    if !isdir(paths_to_tests[i])
        mkdir(paths_to_tests[i])
    end
    cd(paths_to_tests[i]) do
        mpi = false
        grav = false
        cool = false
        if tests_modules[i] == "Gravity" || tests_modules[i] == "Hydro-Gravity"
            grav = true
        end
        if tests_modules[i] == "Cooling"
            cool = true
        end
        _create_makefile(mpi, grav, cool)
        #_create_parameters
        #_create_users
    end
end
