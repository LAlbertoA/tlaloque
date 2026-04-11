import REPL
using REPL.TerminalMenus

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

options = ["Run all tests", "Run tests for a specific module", "Run specific tests", "Exit"]
menu = RadioMenu(options, pagesize = length(options))
choice = request("Select an option:", menu)

test_module_options = ["Hydrodynamics", "Hydrodynamics with Cooling", "Gravity", "Hydrodynamics with Gravity"]

test_hydro_options = ["Sod", "Sedov-Taylor", "Blast Wave", "Kelvin-Helmholtz Instability", "Rayleigh-Taylor Instability", "Orszag-Tang Vortex",
    "Gresho Vortex", "Noh Problem", "Double Mach Reflection", "Wind Tunnel with a Step"]
test_hydro_cooling_options = ["Radiative Sedov-Taylor"]
test_gravity_options = ["Gravitational Potential Accuracy Test (Multigrid)", "Point Potential Test"]
test_hydro_gravity_options = ["Jeans Instability Test", "Gravitational Collapse Test", "Evrard's collapse", "Truelove Collapse Test"]

if choice == 1
    println("Option ", options[choice], " selected!")
    println("Running all tests...")
    println("Current tests in this option:")
    for i in 1:length(test_hydro_options)
        println(" - ", test_hydro_options[i])
    end
    for i in 1:length(test_hydro_cooling_options)
        println(" - ", test_hydro_cooling_options[i])
    end
    for i in 1:length(test_gravity_options)
        println(" - ", test_gravity_options[i])
    end
    for i in 1:length(test_hydro_gravity_options)
        println(" - ", test_hydro_gravity_options[i])
    end
    println("Running ...")
elseif choice == 2
    println("Option ", options[choice], " selected!")
    test_module_menu = RadioMenu(test_module_options, pagesize = length(test_module_options))
    test_module_choice = request("Select a module to test:", test_module_menu)
    if test_module_choice == 1
        println("Running tests for Hydrodynamics...")
        println("Current tests in this module:")
        for i in 1:length(test_hydro_options)
            println(" - ", test_hydro_options[i])
        end
    elseif test_module_choice == 2
        println("Running tests for Hydrodynamics with Cooling...")
        println("Current tests in this module:")
        for i in 1:length(test_hydro_cooling_options)
            println(" - ", test_hydro_cooling_options[i])
        end
    elseif test_module_choice == 3
        println("Running tests for Gravity...")
        println("Current tests in this module:")
        for i in 1:length(test_gravity_options)
            println(" - ", test_gravity_options[i])
        end
    elseif test_module_choice == 4
        println("Running tests for Hydrodynamics with Gravity...")
        println("Current tests in this module:")
        for i in 1:length(test_hydro_gravity_options)
            println(" - ", test_hydro_gravity_options[i])
        end
    else
        println("Invalid choice. Please select a valid module.")
    end
elseif choice == 3
    println("Option ", options[choice], " selected!")
    println("Running specific tests...")
    println("Current tests available:")
    println("Menu canceled.")
end