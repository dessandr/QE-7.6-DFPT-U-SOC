!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
! Copyright (C) 2021-2023 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Ivan Maliyov
!                         Dhruv Desai, Sergio Pineda Flores, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  main entrance for the perturbo code.
!
! Maintenance:
!===============================================================================

program qe2pert
   use mp_global, only: mp_startup, mp_global_end
   use qe_mpi_mod, only: ionode
   use input_param, only: read_input_param, load_ephmat, debug, lwannier, &
      yaml_fname
   use electronic_data, only: init_electronic_data, deallocate_evc_bec
   use lattice_data, only: init_lattice_data
   use elph_matrix, only: compute_elph_matrix, validate_soc_u_load_only, &
      validate_soc_u_production
   use perturbo_data, only: save_perturbo_data
   use yaml_utils,  only: open_output_yaml, close_output_yaml
   use yaml_utils_qe2pert, only: output_param_to_yaml_qe2pert
   use ldaU, only: lda_plus_u
   use noncollin_module, only: noncolin
   implicit none
   logical :: soc_u_load_only
   integer :: load_only_env_status, load_only_env_length
   character(len=256) :: load_only_env

   call mp_startup()
   call environment_setup("QE2PERT", .false.)
      
   call read_input_param()

   load_only_env = ''
   call get_environment_variable('QE2PERT_SOC_U_LOAD_ONLY', &
      load_only_env, length=load_only_env_length, &
      status=load_only_env_status)
   soc_u_load_only = .false.
   if(load_only_env_status == -1 .or. &
      load_only_env_length > len(load_only_env)) &
      call errore('qe2pert', &
         'QE2PERT_SOC_U_LOAD_ONLY value is too long', 1)
   if(load_only_env_status /= 0 .and. load_only_env_status /= 1) &
      call errore('qe2pert', &
         'cannot read QE2PERT_SOC_U_LOAD_ONLY', 1)
   if(load_only_env_status == 0) then
      select case(trim(adjustl(load_only_env)))
      case('1', 'true', 'TRUE', 'yes', 'YES')
         soc_u_load_only = .true.
      case('0', 'false', 'FALSE', 'no', 'NO')
         soc_u_load_only = .false.
      case default
         call errore('qe2pert', &
            'invalid QE2PERT_SOC_U_LOAD_ONLY value', 1)
      end select
   endif
   if(soc_u_load_only .and. &
      ((.not. debug) .or. lwannier .or. load_ephmat)) &
      call errore('qe2pert', 'SOC+U load-only requires debug=.true., ' // &
         'lwannier=.false., and load_ephmat=.false.', 1)

   ! Open the YAML file
   if(ionode) then
      call open_output_yaml(yaml_fname, 'qe2pert')
      call output_param_to_yaml_qe2pert()
   end if
 
   ! read prefix.save from pwscf calculation.
   call init_electronic_data()
   ! read lattice dynamical from PHonon calculations
   call init_lattice_data()

   if(soc_u_load_only .and. ((.not. lda_plus_u) .or. (.not. noncolin))) &
      call errore('qe2pert', &
         'SOC+U load-only requires a noncollinear Hubbard restart', 1)

   if( soc_u_load_only ) then
      call validate_soc_u_load_only()
   else
      if(lda_plus_u .and. noncolin .and. (.not. load_ephmat)) &
         call validate_soc_u_production()
      call save_perturbo_data()

      if( load_ephmat ) then
         call save_elph_mat_wann()
      else
         call compute_elph_matrix()
         !release large arrays that will not be used.
         call deallocate_evc_bec()
         !output data to hdf5 file
         call collect_elph_matrix()
         ! compute and save g(Re, Rp)
         if(.not. debug) call save_elph_mat_wann()
      endif
   endif
   
   call environment_stop("QE2PERT")

   ! Close the YAML file
   if(ionode) then
      call close_output_yaml()
   end if

   call mp_global_end()
end program qe2pert
