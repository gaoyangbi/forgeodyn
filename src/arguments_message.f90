    ! module pygeodyn_message
    ! implicit none
    
    ! contains
    ! print message
    subroutine print_version()
        write (*,'(A)') "**********************************************"
        write (*,'(A)',advance='no') "pygeodyn" // CHAR(10) // &
        "Copyright (C) 2019 Geodynamo"// CHAR(10) // &
            CHAR(10) // &
        "This program is free software: you can redistribute it and/or modify" // CHAR(10) // &
        "it under the terms of the GNU General Public License as published by"// CHAR(10) // &
        "the Free Software Foundation, either version 3 of the License, or"// CHAR(10) // &
        "(at your option) any later version."// CHAR(10) // &
            CHAR(10) // &
        "This program is distributed in the hope that it will be useful,"// CHAR(10) // &
        "but WITHOUT ANY WARRANTY; without even the implied warranty of"// CHAR(10) // &
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"// CHAR(10) // &
        "GNU General Public License for more details."// CHAR(10) // &
            CHAR(10) // &
        "You should have received a copy of the GNU General Public License"// CHAR(10) // &
        "along with this program.  If not, see <http://www.gnu.org/licenses/>."// CHAR(10) // &
            CHAR(10) // &
    
        "This FORTRAN code of the program is modified by YG Zhang in Wuhan University"// CHAR(10) // &
        "3/12/2025"// CHAR(10)
        write (*,'(A)') "Program version: 1.0.0"
        write (*,'(A)') "**********************************************"// CHAR(10)
    end subroutine print_version

    ! parameters message 
    subroutine print_help()
        write (*,'(A)') "Help information:"
        write (*,'(A)') "  -h, --help     Display this help message"
        write (*,'(A)') "  -v, --version  Display the program version"
        write (*,'(A)') "**********************************************"
        write (*,'(A)') " -conf,   type=str, default=con_file,  help: config file"
        write (*,'(A)') " -m,      type=int, default=20,            help: Number of realisations to consider"
        write (*,'(A)') " -algo,   type=str, default='augkf',      help: Algorithm to use"
        write (*,'(A)') " -shear,  type=int, default=0,            help: Compute shear at analysis time: 0 NO, 1 YES"
        write (*,'(A)') " -seed,   type=int, default=0,         help: Seed for initializing random states"
        write (*,'(A)') " -d,      type=int, default=2,            help: Logging level: 1 DEBUG, 2 INFO, 3 WARNING, 4 ERROR, 5 CRITICAL, default=INFO"
        write (*,'(A)') " -path,   type=str, default='./pygeodyn_results/', help: Path of the output folder"
        write (*,'(A)') " -cname,  type=str, default='Current_computation', help: Name of the computation"
        write (*,'(A)') " -l,      type=str, default='log',                 help: Name of the log file"
        write (*,'(A)') "**********************************************"// CHAR(10)
    end subroutine print_help
    
    subroutine print_usage()
        write (*,'(A)') "Usage: geodyn_fortran [-h | --help] [-v | --version] [-conf ......]"// CHAR(10)
    end subroutine print_usage

    
    ! end module