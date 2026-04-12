function _run_menu()

    ## Function that creates the main menu for the tests
    ## This is the first Command line that i create so probably very weird but it works
    ## First options to see which tests wants the user to run
    options = ["Run all tests", "Run tests for a specific module", "Run specific tests", "Exit"]
    menu = RadioMenu(options, pagesize = length(options))
    choice = request("Select an option:", menu)

    ## Module options
    test_module_options = ["Hydrodynamics", "Hydrodynamics with Cooling", "Gravity", "Hydrodynamics with Gravity"]

    ## Tests in each module
    test_hydro_options = ["Sod", "Sedov-Taylor", "Blast Wave", "Kelvin-Helmholtz Instability", "Rayleigh-Taylor Instability", "Orszag-Tang Vortex",
        "Gresho Vortex", "Noh Problem", "Double Mach Reflection", "Wind Tunnel with a Step"]
    test_hydro_cooling_options = ["Radiative Sedov-Taylor"]
    test_gravity_options = ["Gravitational Potential Accuracy Test (Multigrid)", "Point Potential Test"]
    test_hydro_gravity_options = ["Jeans Instability Test", "Gravitational Collapse Test", "Evrard's collapse", "Truelove Collapse Test"]

    ## All tests
    tests = cat("Hydro: " .* test_hydro_options, 
                "Cooling: " .* test_hydro_cooling_options,
                "Gravity: " .* test_gravity_options,
                "Hydro-Gravity: " .* test_hydro_gravity_options,
                dims=1)

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
        return tests, trues(length(tests))
    elseif choice == 2
        println("Option ", options[choice], " selected!")
        test_module_menu = RadioMenu(test_module_options, pagesize = min(length(test_module_options),8))
        test_module_choice = request("Select a module to test:", test_module_menu)
        if test_module_choice == 1
            println("Running tests for Hydrodynamics...")
            println("Current tests in this module:")
            for i in 1:length(test_hydro_options)
                println(" - ", test_hydro_options[i])
            end
            return tests, [startswith(i, "Hydro:") for i in tests]
        elseif test_module_choice == 2
            println("Running tests for Hydrodynamics with Cooling...")
            println("Current tests in this module:")
            for i in 1:length(test_hydro_cooling_options)
                println(" - ", test_hydro_cooling_options[i])
            end
            return tests, [startswith(i, "Cooling:") for i in tests]
        elseif test_module_choice == 3
            println("Running tests for Gravity...")
            println("Current tests in this module:")
            for i in 1:length(test_gravity_options)
                println(" - ", test_gravity_options[i])
            end
            return tests, [startswith(i, "Gravity:") for i in tests]
        elseif test_module_choice == 4
            println("Running tests for Hydrodynamics with Gravity...")
            println("Current tests in this module:")
            for i in 1:length(test_hydro_gravity_options)
                println(" - ", test_hydro_gravity_options[i])
            end
            return tests, [startswith(i, "Hydro-Gravity:") for i in tests]
        else
            println("Invalid choice. Please select a valid module.")
        end
    elseif choice == 3
        println("Option ", options[choice], " selected!")
        println("Select specific tests to run:")
        tests_menu = MultiSelectMenu(tests, pagesize = min(length(tests),8))
        tests_choice = request("Select tests to run:", tests_menu)

        println("Tests selected: ")
        for test in tests_choice
            println(" - ", tests[test])
        end
        return tests, [i in tests_choice for i in 1:length(tests)]
    end
end

function _create_makefile(mpi,gravity,cooling)
    
    mpi_string = "N"
    gravity_string = "N"
    cooling_string = "N"
    
    if mpi==true
        mpi_string = "Y"
    end
    if gravity==true
        gravity_string = "Y"
    end
    if cooling==true
        cooling_string = "Y"
    end
    open("Makefile", "w") do f
        write(f, """
            PROGRAM= tlaloque
            COMPILER= gfortran
            ######  USER_FLAGS for gfortran:
            #USER_FLAGS= -O3 -Wall -fcheck=all
            ######  USER_FLAGS for ifort:
            USER_FLAGS= -O3 #-warn all -check all
            MPI= $mpi_string
            DOUBLEP= Y
            PASB= N
            GRV= $gravity_string
            COOL= $cooling_string

            MODULES_USER= \\

            MODULES_MAIN= \\
            ./source/constants.o \\
            ./user/parameters.o \\
            ./source/globals.o \\
            ./source/cooling.o \\
            ./source/pointsource_gravity.o \\
            ./source/sources.o \\
            ./source/winds.o \\
            ./user/user.o 

            OBJECTS_MAIN= \\
            ./source/HydroCore.o \\
            ./source/init.o \\
            ./source/boundaries.o \\
            ./source/HLLC.o \\
            ./source/mainHD.o 

            CFLAGS = \$(USER_FLAGS) -cpp

            ifeq (\$(MPI),Y)
            CFLAGS += -DMPIP
            endif

            ifeq (\$(COOL),Y)
            CFLAGS += -DCOOL
            endif

            ifeq (\$(GRV),Y)
            CFLAGS += -DGRAV
            endif

            ifeq (\$(DOUBLEP),Y)
            CFLAGS += -DDOUBLEP
            ifeq (\$(COMPILER),ifort)
            CFLAGS += -r8
            endif
            ifeq (\$(COMPILER),gfortran)
            CFLAGS += -fdefault-real-8
            endif
            endif

            ifeq (\$(PASB),Y)
            CFLAGS += -DPASBP
            endif

            ifeq (\$(MPI),Y)
            COMPILER = mpif90
            endif

            OBJECTS_ALL = \${MODULES_MAIN} \${MODULES_USER} \${OBJECTS_MAIN}

            \$(PROGRAM) : prebuild \${OBJECTS_ALL}
                @echo Linking object files ...
                @\$(COMPILER) \$(CFLAGS) \${OBJECTS_ALL} -o \$@
                @echo Cleaning up ...
                @rm -f *.o *.mod source/*.o source/*.mod user/*.o user/*.mod
                @echo "Done! (`date`)"

            prebuild :
                @echo "\$(PROGRAM) build started `date`"

            %.o : %.f90
                @echo Compiling \$^ ...
                @\$(COMPILER) \$(CFLAGS) -c \$^ -o \$@

            clean :
                rm -f *.o *.mod source/*.o source/*.mod user/*.o user/*.mod

            cleanall :
                rm -f *.o *.mod source/*.o source/*.mod user/*.o user/*.mod
                rm -f \$(PROGRAM).*
                rm -f coldens
                rm -f extract
                rm -f DATA/*.bin
                rm -f DATA/*.vtk
                rm -f DATA/*.dat
                rm -f DATA/*.log
                rm -f DATA/*.visit"""
        )
    end
end
