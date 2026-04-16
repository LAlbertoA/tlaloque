import REPL
using REPL.TerminalMenus

include("tests_utilities.jl")

main_path = @__DIR__
tests_path = joinpath(main_path,"tests")

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
println("Running this tests will make a new folder in the tlaloque directory called 'tests' and one
folder inside it for each test. Each folder test will be named after the test. Inside each individual
test folder there will be a Makefile, the user.f90 and parameters.f90 files of the test, and a DATA
folder with the outputs of the test. Also, the test will be compared with expected results and this will
be saved to a file called '(test name)_results.txt'.\n")

println("THIS PROCESS CAN TAKE A WHILE, DEPENDING ON THE NUMBER OF TESTS YOU SELECTED AND YOUR COMPUTER PERFORMANCE.
IF YOU WANT TO CONTINUE, SELECT 'Continue with tests' IN THE NEXT MENU. OTHERWISE, PRESS CTRL + C
TO ABORT, SELECT 'Abort' OR PRESS Q IN THE NEXT MENU.")
println("\n ################################################################################### \n")

continue_options = ["Continue with tests", "Abort"]
continue_options_menu = RadioMenu(continue_options, pagesize = 2)
continue_options_selected = request("Select an option:", continue_options_menu)
if continue_options_selected == 2
    println("Tests Aborted.")
    exit()
elseif continue_options_selected == -1
    println("Tests Aborted.")
    exit()
end

if !isdir(tests_path)
    mkdir(tests_path)
end

println("Continuing with tests. `tests` folder created.\n")

println("This process can take a while, depending on the number of tests you selected and your computer performance.\n")

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
        tname = "tlaloque_"*split(paths_to_tests[i], "/")[end]
        if tests_modules[i] == "Gravity" || tests_modules[i] == "Hydro-Gravity"
            grav = true
        end
        if tests_modules[i] == "Cooling"
            cool = true
        end
        println("Running test: ", tests_names[i])
        println("####################################################################################")
        println("Creating Makefile, parameters.f90 and user.f90 files in ", paths_to_tests[i], "\n")
        _create_makefile(mpi, grav, cool, paths_to_tests[i], tname)
        _create_parameters(tests_names[i])
        _create_user(tests_names[i])
        println("Files created. Compiling and running the test.\n")
        run(`make clean`)
        run(`make`)
        println("Checking for DATA folder and creating it if it doesn't exist.\n")
        if !isdir("DATA")
            mkdir("DATA")
        end
        println("Running the test and saving outputs to DATA folder.\n")
        run(pipeline(`./$tname`, "DATA/output.txt"))
        println("Test run completed. Comparing results with expected results.\n")
    end
end
