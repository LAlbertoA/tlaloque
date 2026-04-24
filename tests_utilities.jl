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
    test_hydro_options = ["Sod", "Sedov-Taylor"]
    ## Future hydro tests:  "Blast Wave", "Kelvin-Helmholtz Instability", "Rayleigh-Taylor Instability", "Orszag-Tang Vortex",
    ##    "Gresho Vortex", "Noh Problem", "Double Mach Reflection", "Wind Tunnel with a Step"
    test_hydro_cooling_options = ["Radiative Sedov-Taylor"]
    test_gravity_options = ["Gravitational Potential Accuracy Test (Multigrid)", "Point Potential Test"]
    test_hydro_gravity_options = ["Evrard's Collapse", "Truelove Collapse Test"]
    ## Future hydro-gravity tests: "Jeans Instability Test", "Gravitational Collapse Test",
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
    elseif choice == 4
        println("Exiting.")
        exit()
    elseif choice == -1
        println("Exiting.")
        exit()
    end
end

function _create_makefile(mpi,gravity,cooling,tpath,tname)
    
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
            PROGRAM= $tname
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
            $main_path/source/constants.o \\
            $tpath/parameters.o \\
            $main_path/source/globals.o \\
            $main_path/source/cooling.o \\
            $main_path/source/pointsource_gravity.o \\
            $main_path/source/sources.o \\
            $main_path/source/winds.o \\
            $tpath/user.o 

            OBJECTS_MAIN= \\
            $main_path/source/HydroCore.o \\
            $main_path/source/init.o \\
            $main_path/source/boundaries.o \\
            $main_path/source/HLLC.o \\
            $main_path/source/mainHD.o 

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
            \t@echo Linking object files ...
            \t@\$(COMPILER) \$(CFLAGS) \${OBJECTS_ALL} -o \$@
            \t@echo Cleaning up ...
            \t@rm -f $main_path/*.o $main_path/*.mod $main_path/source/*.o $main_path/source/*.mod $tpath/*.o $tpath/*.mod
            \t@echo "Done! (`date`)"

            prebuild :
            \t@echo "\$(PROGRAM) build started `date`"

            %.o : %.f90
            \t@echo Compiling \$^ ...
            \t@\$(COMPILER) \$(CFLAGS) -c \$^ -o \$@

            clean :
            \trm -f $main_path/*.o $main_path/*.mod $main_path/source/*.o $main_path/source/*.mod $tpath/*.o $tpath/*.mod

            cleanall :
            \trm -f $main_path/*.o $main_path/*.mod $main_path/source/*.o $main_path/source/*.mod $tpath/*.o $tpath/*.mod
            \trm -f \$(PROGRAM).*
            \trm -f coldens
            \trm -f extract
            \trm -f DATA/*.bin
            \trm -f DATA/*.vtk
            \trm -f DATA/*.dat
            \trm -f DATA/*.log
            \trm -f DATA/*.visit"""
        )
    end
end

function _parameters(name)

    if name == "Sod"
        boundaries = ["OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"]
        cells = [150,150,150]
        physical_size = [-0.25, 0.25, -0.25, 0.25, -0.25, 0.25]
        gamma = "5.0/3.0"
        nout = 2
        size = 0
        Gconst = "GR"
        cfl = 0.8
        eta = 0.9e-1
        tfin = "0.12"
    elseif name == "Sedov-Taylor"
        boundaries = ["OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"]
        cells = [150,150,150]
        physical_size = [-5.0, 5.0, -5.0, 5.0, -5.0, 5.0]
        gamma = "5.0/3.0"
        nout = 2
        size = 0
        Gconst = "GR"
        cfl = 0.8
        eta = 0.9e-1
        tfin = "0.1"
    elseif name == "Radiative Sedov-Taylor"
        boundaries = ["OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"]
        cells = [150,150,150]
        physical_size = ["-30.0*PC", "30.0*PC", "-30.0*PC", "30.0*PC", "-30.0*PC", "30.0*PC"]
        gamma = "5.0/3.0"
        nout = 10
        size = 0
        Gconst = "GR"
        cfl = 0.8
        eta = 0.9e-1
        tfin = "50*YEAR*1.e3"
    elseif name == "Gravitational Potential Accuracy Test (Multigrid)"
        boundaries = ["OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"]
        cells = [100,100,100]
        physical_size = [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5]
        gamma = "5.0/3.0"
        nout = 1
        size = 0
        Gconst = "GR"
        cfl = 0.8
        eta = 0.5e-2
        tfin = "0"
    elseif name == "Evrard's Collapse"
        boundaries = ["REFLECTIVE", "REFLECTIVE", "REFLECTIVE", "REFLECTIVE", "REFLECTIVE", "REFLECTIVE"]
        cells = [100,100,100]
        physical_size = [-2.0, 2.0, -2.0, 2.0, -2.0, 2.0]
        gamma = "5.0/3.0"
        nout = 10
        size = 0
        Gconst = 1.0
        cfl = 0.8
        eta = 0.5e-2
        tfin = "3"
    elseif name == "Truelove Collapse Test"
        boundaries = ["OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"]
        cells = [100,100,100]
        physical_size = ["-535.0*AU", "535.0*AU", "-535.0*AU", "535.0*AU", "-535.0*AU", "535.0*AU"]
        gamma = "1.001"
        nout = 10
        size = 0
        Gconst = "GR"
        cfl = 0.8
        eta = 0.5e-2
        tfin = "2110*YEAR"
    end

    return boundaries, cells, physical_size, gamma, nout, size, Gconst, cfl, eta, tfin

end

function _create_parameters(name)

    boundaries, cells, physical_size, gamma, nout, size, Gconst, cfl, eta, tfin = _parameters(name)

    open("parameters.f90", "w") do f
        write(f, """
        module parameters

            use constants
            implicit none

        #ifdef MPIP
            include "mpif.h"
        #endif
                                                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            integer, parameter  :: LB = BC_$(boundaries[1])    !!                           Boundaries                           !!
            integer, parameter  :: RB = BC_$(boundaries[2])    !!  LB = Left Boundary   DB = Down Boundary  FB = Front Boundary  !!
            integer, parameter  :: TB = BC_$(boundaries[3])    !!  RB = Right Boundary  TB = Top Boundary   BB = Back Boundary   !!
            integer, parameter  :: DB = BC_$(boundaries[4])    !!           (x)                 (y)                  (z)         !!
            integer, parameter  :: FB = BC_$(boundaries[5])    !!                        Type of boundary                        !!
            integer, parameter  :: BB = BC_$(boundaries[6])    !!  1 = Outflow    2 = Reflective    3 = Periodic    4 = Inflow   !!
                                                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            integer, parameter  :: limtr = LIMITER_VANLEER  !! Limiter

            integer, parameter  :: ndim = 3           !! Number of dimensions of the problem (Currently only supporting 3).
            integer, parameter  :: neq = 5            !! Number of equations to solve (Currently only supporting 5).
            integer, parameter  :: size = $size           !! 
            integer, parameter  :: nxtot = $(cells[1])  !! Number of cells in each axis. If using self-gravity, nx, ny and nz
            integer, parameter  :: nytot = $(cells[2])  !! MUST be a power of 2 with `size` it's power: nx, ny, nz = 2**size.
            integer, parameter  :: nztot = $(cells[3])  !! If using MPI and self-gravity, nx, ny and nz MUST be a power of 2
            integer, parameter  :: choice = 2         !! like so: nx, ny, nz = (mpix, mpiy, mpiz)*2**size

            integer, parameter      :: num_out = $nout  !! Number of output files 
            integer, parameter      :: outfile = VTK  !! Type of output. Options are VTK, DAT, BIN
            character(*), parameter :: outputpath = './DATA/'
        #ifdef COOL
            character(*), parameter :: cooling_file = '../../Z1.0.dat' !! Cooling table file location for atomic cooling
        #endif
            character(*), parameter :: posfile = './posest75.dat' !! Star positions file for winds and point_gravity modules
            integer, parameter      :: nghost = 2   !! Order

            real, parameter  :: xl = $(physical_size[1])           !! Position of first physical cell in the x axis
            real, parameter  :: xr = $(physical_size[2])           !! Position of last physical cell in the x axis
            real, parameter  :: yl = $(physical_size[3])           !! Position of first physical cell in the y axis
            real, parameter  :: yr = $(physical_size[4])           !! Position of last physical cell in the y axis
            real, parameter  :: zl = $(physical_size[5])           !! Position of first physical cell in the z axis
            real, parameter  :: zr = $(physical_size[6])           !! Position of last physical cell in the z axis
            real, parameter  :: gamma = $gamma        !! Specific heat ratio. C_p/C_v
            real, parameter  :: mu0 = 1.0              !! Mean particule mass
            real, parameter  :: Gconst = $Gconst       !! Gravitational constant. If using physical units, GR = 6.67430e-8 dyn cm^2/g^2
            real, parameter  :: cfl = $cfl              !! Courant-Friedrichs-Lewy condition. 0 < cfl < 1
            real, parameter  :: eta = $eta              !! Artificial Viscosity. 0 < eta < 0.5
            real, parameter  :: tfin = $tfin            !! Total simulated time.
            real, parameter  :: cooling_limit = 0.5
                        !!! MPI Constants definitions !!!
        #ifdef MPIP          
            integer, parameter  :: mpix = 4            !! Number of subdivisions of the computational domain in the x axis
            integer, parameter  :: mpiy = 4            !! Number of subdivisions of the computational domain in the y axis
            integer, parameter  :: mpiz = 4            !! Number of subdivisions of the computational domain in the z axis 
            integer, parameter  :: np = mpix*mpiy*mpiz !! Total number of processor cores needed for the given subdivisions
            integer, parameter  :: nx = nxtot/mpix     !! Number of cells in each core in the x axis
            integer, parameter  :: ny = nytot/mpiy     !! Number of cells in each core in the y axis
            integer, parameter  :: nz = nztot/mpiz     !! Number of cells in each core in the z axis
        #else
            integer, parameter  :: mpix = 1
            integer, parameter  :: mpiy = 1
            integer, parameter  :: mpiz = 1
            integer, parameter  :: nx = nxtot
            integer, parameter  :: ny = nytot
            integer, parameter  :: nz = nztot
            integer, parameter  :: np = 1
        #endif

        #ifdef MPIP
        #ifdef DOUBLEP
            integer, parameter :: mpi_real_kind=mpi_real8
        #else
            integer, parameter :: mpi_real_kind=mpi_real4
        #endif
        #else
            integer, parameter :: mpi_real_kind=0
        #endif

            real, parameter  :: dtout = tfin/num_out   !! Time step between outputs

            integer, parameter  :: nxmax = nx+nghost    !!
            integer, parameter  :: nymax = ny+nghost    !!
            integer, parameter  :: nzmax = nz+nghost    !! Total number of cells including
                                                        !! ghost cells in each core given
            integer, parameter  :: nxmin = 1-nghost     !! the order and reconstruction scheme
            integer, parameter  :: nymin = 1-nghost     !!
            integer, parameter  :: nzmin = 1-nghost     !!

            real, parameter     :: dx = (xr-xl)/nxtot   !! cell size in the x axis
            real, parameter     :: dy = (yr-yl)/nytot   !! cell size in the y axis
            real, parameter     :: dz = (zr-zl)/nztot   !! cell size in the z axis
        #ifdef GRAV
        #ifdef MPIP
            integer, parameter  :: lvlm = size
        #else
            integer, parameter  :: lvlm = size
        #endif
        #endif
            logical, parameter      :: logged = .false.            !! Sets whether or not create a logging file
            integer, parameter      :: logu = 12
            character(*), parameter :: logfile='./log125SGG.dat'   !! Name of the logging file

            character(*), parameter :: paramFile = "tlaloque.dat"  !! Name of the file to print all runtime parameters
        end module parameters
        """)
    end
end

function _write_initconds(name)

    if name == "Sod"
        initconds_name = "initflowSod"
        init_conds = """
        subroutine initflowSod()

            use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma
            use globals
            use constants

            implicit none

            real, parameter       :: rhoO = 4.0, uO = 0.0, vO = 0.0, wO = 0.0, pO = 1.0
            real, parameter       :: rhoI = 1.0, uI = 0.0, vI = 0.0, wI = 0.0, pIn = 0.1795
            real, parameter       :: xc = 0.0, yc = 0.0, zc = 0.0
            integer               :: i, j, k
            real                  :: x, y, z

            x = 0.0
            y = 0.0
            z = 0.0

            do i = nxmin, nxmax
                do j = nymin, nymax
                    do k = nzmin, nzmax
                        call xcoord(i, x)
                        call xcoord(j, y)
                        call xcoord(k, z)
                        if (-x-y-z > 0.0) then
                            U(1,i,j,k) = rhoO
                            U(2,i,j,k) = rhoO*uO
                            U(3,i,j,k) = rhoO*vO
                            U(4,i,j,k) = rhoO*wO
                            U(5,i,j,k) = 0.5*rhoO*(uO**2 + vO**2 + wO**2)+pO/(gamma-1.0)
                        else
                            U(1,i,j,k) = rhoI
                            U(2,i,j,k) = rhoI*uI
                            U(3,i,j,k) = rhoI*vI
                            U(4,i,j,k) = rhoI*wI
                            U(5,i,j,k) = 0.5*rhoI*(uI**2 + vI**2 + wI**2)+pIn/(gamma-1.0)
                        endif
                    enddo
                enddo
            enddo
        #ifdef GRAV
            PHI(:,:,:) = 0.0
        #endif
        end subroutine initflowSod
        """
    elseif name == "Sedov-Taylor"
        initconds_name = "initflowSedov"
        init_conds = """
        subroutine initflowSedov()

            use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma, dx
            use globals
            use constants

            implicit none

            real, parameter       :: rhoO = 1.0, uO = 0.0, vO = 0.0, wO = 0.0, pO = 1.5e-20*(gamma-1.0)
            real, parameter       :: rhoI = 1.0, uI = 0.0, vI = 0.0, wI = 0.0
            real, parameter       :: xc = 0.0, yc = 0.0, zc = 0.0, pIn = 2.481617647e6*(gamma-1.0)
            integer               :: i, j, k
            real                  :: x, y, z, r

            x = 0.0
            y = 0.0
            z = 0.0

            do i = nxmin, nxmax
                do j = nymin, nymax
                    do k = nzmin, nzmax
                        call xcoord(i, x)
                        call ycoord(j, y)
                        call zcoord(k, z)
                        r = sqrt(x**2+y**2+z**2)
                        if (r < 0.2) then
                            U(1,i,j,k) = rhoI
                            U(2,i,j,k) = rhoI*uI
                            U(3,i,j,k) = rhoI*vI
                            U(4,i,j,k) = rhoI*wI
                            U(5,i,j,k) = 0.5*rhoI*(uI**2 + vI**2 + wI**2)+pIn/(gamma-1.0)
                        else
                            U(1,i,j,k) = rhoO
                            U(2,i,j,k) = rhoO*uO
                            U(3,i,j,k) = rhoO*vO
                            U(4,i,j,k) = rhoO*wO
                            U(5,i,j,k) = 0.5*rhoO*(uO**2 + vO**2 + wO**2)+pO/(gamma-1.0)
                        endif
                    enddo
                enddo
            enddo
        #ifdef GRAV
            PHI(:,:,:) = 0.0
        #endif
        end subroutine initflowSedov
        """
    elseif name == "Radiative Sedov-Taylor"
        initconds_name = "initflowRadiativeSedov"
        init_conds = """
        subroutine initflowRadiativeSedov()

            use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma
            use globals
            use constants

            implicit none

            real, parameter       :: rhoO = AMU*1.0, uO = 0.0, vO = 0.0, wO = 0.0, pO = (1.0*AMU*KB*10**4)/AMU
            real, parameter       :: rhoI = AMU*1.0, uI = 0.0, vI = 0.0, wI = 0.0, E0 = 0.5e51, r0 = 0.35*PC
            real, parameter       :: xc = 0.0, yc = 0.0, zc = 0.0, pIn = (E0/((4.0/3.0)*pi*r0**3))*(gamma-1.0)
            integer               :: i, j, k
            real                  :: x, y, z, r

            x = 0.0
            y = 0.0
            z = 0.0

            do i = nxmin, nxmax
                do j = nymin, nymax
                    do k = nzmin, nzmax
                        call xcoord(i, x)
                        call ycoord(j, y)
                        call zcoord(k, z)
                        r = sqrt(x**2+y**2+z**2)
                        if (r > r0) then
                            U(1,i,j,k) = rhoO
                            U(2,i,j,k) = rhoO*uO
                            U(3,i,j,k) = rhoO*vO
                            U(4,i,j,k) = rhoO*wO
                            U(5,i,j,k) = 0.5*rhoO*(uO**2 + vO**2 + wO**2)+pO/(gamma-1.0)
                        else
                            U(1,i,j,k) = rhoI
                            U(2,i,j,k) = rhoI*uI
                            U(3,i,j,k) = rhoI*vI
                            U(4,i,j,k) = rhoI*wI
                            U(5,i,j,k) = 0.5*rhoI*(uI**2 + vI**2 + wI**2)+pIn/(gamma-1.0)
                        endif
                    enddo
                enddo
            enddo
        #ifdef GRAV
            PHI(:,:,:) = 0.0
        #endif
        end subroutine initflowRadiativeSedov
        """
    elseif name == "Gravitational Potential Accuracy Test (Multigrid)"
        initconds_name = "initflowGravPotAccTest"
        init_conds = """
        subroutine initflowGravPotAccTest()

            use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma
            use globals
            use constants

            implicit none

            real, parameter       :: rhoI = 1.0, uI = 0.0, vI = 0.0, wI = 0.0, pI = 1.5e-5*(gamma-1.0)
            real, parameter       :: rho0 = 0.0, u0 = 0.0, v0 = 0.0, w0 = 0.0, p0 = 1.5e-5*(gamma-1.0)
            real, parameter       :: xc = 0.0, yc = 0.0, zc = 0.0, r0 = 0.25
            integer               :: i, j, k
            real                  :: x, y, z, r

            x = 0.0
            y = 0.0
            z = 0.0

            do i = nxmin, nxmax
                do j = nymin, nymax
                    do k = nzmin, nzmax
                        call xcoord(i, x)
                        call ycoord(j, y)
                        call zcoord(k, z)
                        r = sqrt(x**2+y**2+z**2)
                        if (r > r0) then
                            U(1,i,j,k) = rho0
                            U(2,i,j,k) = rho0*u0
                            U(3,i,j,k) = rho0*v0
                            U(4,i,j,k) = rho0*w0
                            U(5,i,j,k) = 0.5*rho0*(u0**2 + v0**2 + w0**2)+p0/(gamma-1.0)
                        else
                            U(1,i,j,k) = rhoI
                            U(2,i,j,k) = rhoI*uI
                            U(3,i,j,k) = rhoI*vI
                            U(4,i,j,k) = rhoI*wI
                            U(5,i,j,k) = 0.5*rhoI*(uI**2 + vI**2 + wI**2)+pI/(gamma-1.0)
                        endif
                    enddo
                enddo
            enddo

        #ifdef GRAV
            call pointmass_potential(PHI)
        #endif
        """
    elseif name == "Evrard's Collapse"
        initconds_name = "initflowEvrard"
        init_conds = """
        subroutine initflowEvrard()

            use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma
            use globals
            use constants

            implicit none

            real, parameter       :: uO = 0.0, vO = 0.0, wO = 0.0, pO = (gamma-1.0)*0.05
            real, parameter       :: uI = 0.0, vI = 0.0, wI = 0.0
            integer               :: i, j, k
            real                  :: x, y, z, r, M, R1

            x = 0.0
            y = 0.0
            z = 0.0

            M = 1.0
            R1 = 1.0

            do i = nxmin, nxmax
                do j = nymin, nymax
                    do k = nzmin, nzmax
                        call xcoord(i, x)
                        call ycoord(j, y)
                        call zcoord(k, z)
                        r = sqrt(x**2+y**2+z**2)
                        if (r>R1) then
                            U(1,i,j,k) = 1.e-4
                            U(2,i,j,k) = 0.0
                            U(3,i,j,k) = 0.0
                            U(4,i,j,k) = 0.0
                            U(5,i,j,k) = 1.e-4*0.05
                        else
                            U(1,i,j,k) = M/(2.0*PI*r*R1**2)
                            U(2,i,j,k) = 0.0
                            U(3,i,j,k) = 0.0
                            U(4,i,j,k) = 0.0
                            U(5,i,j,k) = 0.05/(2.0*PI*r*R1**2)
                        end if
                    enddo
                enddo
            enddo
        #ifdef GRAV
            PHI(:,:,:) = 0.0
        #endif
        end subroutine initflowEvrard
        """
    elseif name == "Truelove Collapse Test"
        initconds_name = "initflowTruelove"
        init_conds = """
        subroutine initflowTruelove()

            use parameters, only: nxmax, nymax, nzmax, nxmin, nymin, nzmin, neq, gamma, xr, yr, zr, dx
            use globals
            use constants

            implicit none

            real, parameter       :: xc = 0, yc = 0, zc = 0, cs = 0.167*KMS
            real, parameter       :: uO = 0.0, vO = 0.0, wO = 0.0, rhoO = 0.01*1.0e-15, pO = (rhoO*cs**2)
            real, parameter       :: uI = 0.0, vI = 0.0, wI = 0.0, rhoI = 1.0e-15, pIn = pO
            integer               :: i, j, k
            real                  :: x, y, z, r, R1

            x = 0.0
            y = 0.0
            z = 0.0

            R1 = 7.8e15

            do i = nxmin, nxmax
                do j = nymin, nymax
                    do k = nzmin, nzmax
                        call xcoord(i, x)
                        call ycoord(j, y)
                        call zcoord(k, z)
                        r = sqrt(x**2+y**2+z**2)
                        if (r>R1) then
                            U(1,i,j,k) = rhoO
                            U(2,i,j,k) = 0.0
                            U(3,i,j,k) = 0.0
                            U(4,i,j,k) = 0.0
                            U(5,i,j,k) = pO/(gamma-1)
                        else
                            U(1,i,j,k) = rhoI
                            U(2,i,j,k) = 0.0
                            U(3,i,j,k) = 0.0
                            U(4,i,j,k) = 0.0
                            U(5,i,j,k) = pIn/(gamma-1)
                        end if
                    enddo
                enddo
            enddo
        #ifdef GRAV
            PHI(:,:,:) = 0.0
        #endif
        end subroutine initflowTruelove
        """
    end

    return initconds_name, init_conds
end

function _create_user(name)

    ic_function_name, initconds_function = _write_initconds(name)

    open("user.f90", "w") do f
        write(f, """
        module user

            use parameters
            use globals, only: posstars, Nstrs
        #ifdef GRAV
            use globals, only: PHIP
        #endif
            use winds
            use pgravity

            implicit none

            type(spherical_wind_type), dimension(:), allocatable :: spherical_wind

            contains

            !! Subroutine to set the initial conditions of the problem to be modeled.
            !! If using the winds module, set the parameters and make the call to the
            !! `imposeSphericalWind` subroutine of each star. The positions of stars
            !! read from the `posfile` are stored in the array `posstars` already in
            !! scope of this module
            subroutine user_initconds(U)

                use constants
                use parameters, only: dx
                
                implicit none

                real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(inout) :: U
                integer                                                                    :: i

                call $ic_function_name()

        #ifdef GRAV
                call pointmass_potential(PHIP)
        #endif     
            end subroutine user_initconds

            !! Subroutine to Set the boundary conditions of the problem.
            !! If using the winds module, place here the call to `imposeSphericalWind`
            !! subroutine to each star.
            subroutine user_boundaries(U)

                implicit none

                real, dimension(neq, nxmin:nxmax, nymin:nymax, nzmin:nzmax), intent(inout) :: U
                integer                                                                    :: i

                !!do i = 0,Nstrs-1

                !! call imposeSphericalWind(spherical_wind(i), U)

                !!enddo

            end subroutine user_boundaries

            !! Subroutine to load star positions for the winds and pointsource_gravity module
            !! First line of file should be the number of stars. From second line to last line
            !! should be the positions of the stars, only one star per line and ordered like: `x  y  z`
            !! with blank space between coordinates
            subroutine load_objects()

                implicit none

                integer                           :: i

                open(unit=33,file=posfile,status="old",action="read")

                read(33,*) Nstrs
                
                allocate(spherical_wind(0:Nstrs-1))
                allocate(posstars(0:2,0:Nstrs-1))
                
                do i = 0,Nstrs-1
                read(33,*) posstars(:,i)
                enddo

                close(33)

            end subroutine load_objects

            $initconds_function
        end module user
        """)
    end
end
