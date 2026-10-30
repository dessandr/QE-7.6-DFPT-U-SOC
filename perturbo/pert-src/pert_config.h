! Defining SCAT_FWD causes perturbo to generate a data structure that
! supports forward traversal, i.e. "source element ordered" calculation.
! #define SCAT_FWD

! Defining SCAT_REV causes perturbo to generate a data structure that
! supports reverse traversal, i.e. "target element ordered" calculation.
#define SCAT_REV

! It is not a problem for both SCAT_FWD and SCAT_REV to be defined.
! However, at least one of them needs to be defined.
#if !defined(SCAT_FWD) && !defined(SCAT_REV)
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif


! Defining STORE_G2DT causes the scatter_channels array to store the two
! -g2*dt1 and -g2*dt2 terms, instead of just holding g2.  This increases
! memory requirements of perturbo, but can result in faster computations.
#define STORE_G2DT


! Defining CDYN_USE_ARRAY_REDUCE causes the OpenMP code to use an
! array-reduction instead of "omp atomic" to perform updates to the
! epcol array.  This means that each OpenMP thread has its own copy
! of epcol, which uses more memory, and combining them at the end
! will take time as well.  However, array-reduction is usually much
! faster than "omp atomic".
#define CDYN_USE_ARRAY_REDUCE

! Defining CDYN_SORT_SCAT_TGTS causes the target-oriented code to
! order the scatter-target array elements by decreasing sub-array
! lengths.  For the GPU this seems to greatly improve performance,
! since all threads in a warp will finish working around the same
! time as each other.
#define CDYN_SORT_SCAT_TGTS

! Whether to split different E or T field components into different 
! calls of OpenACC subroutine. In this way, the tensor of 
! scatchan_mkq2nk_col(2,ns,ncomp) will be reduced into (2,ns). Thus, 
! the access to elements (in global memory) will be more efficient. 
! However, the lunch of multiple OpenAcc will have additional overload 
! as a cost.
! #define SPLIT_COMP

! Calculate fermi and bose distribution in advance. The same (n,k) 
! or (m,kq) or (v,q) will get involved in different scatter channels 
! multiple times. The repeatedly calculation of their distribution could 
! be a waste of run time. The memory overhead of this strategy is negligible.
! #define CAL_FERMI_BOSE

! Whether to split nk2mkq and mkq2nk into two separate arrays. Different
! from dynamics, the scattering rate of nk->mkq is different from mkq->nk.
! Splitting both in separate matrices could minimize the time to access the 
! elements in huge global memory of scatchan_mkq2nk_col(2,ns) or 
! scatchan_mkq2nk_col(2,ns,ncomp)
! #define SPLIT_NK_MKQ
