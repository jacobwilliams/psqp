!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/22/2015
!
!  Numeric kind definitions.

    module psqp_kind_module

    use, intrinsic :: iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter :: wp = real32   !! Real working precision [4 bytes]
#elif REAL64
    integer,parameter :: wp = real64   !! Real working precision [8 bytes]
#elif REAL128
    integer,parameter :: wp = real128  !! Real working precision [16 bytes]
#else
    integer,parameter :: wp = real64   !! Real working precision if not specified [8 bytes]
#endif

    integer,parameter,public :: psqp_wp = wp   !! Working real precision

    end module psqp_kind_module
!*****************************************************************************************
