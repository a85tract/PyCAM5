!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_taumol.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.7 $
!     created:   $Date: 2009/10/20 15:08:37 $
!
      module rrtmg_lw_taumol

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc
!      use parkind, only : im => kind_im, rb => kind_r8
      use parrrtm, only : mg, nbndlw, maxxsec, ngptlw
      use rrlw_con, only: oneminus
      use rrlw_wvn, only: nspa, nspb
      use rrlw_vsn, only: hvrtau, hnamtau

      implicit none

      logical :: use_native_taugb10_11_14_lw_impl = .false.
      logical :: taugb10_11_14_lw_impl_selected = .false.
      logical :: taugb10_11_14_lw_entered_logged(3) = .false.
      logical :: use_native_taugb1_2_6_8_12_16_lw_impl = .false.
      logical :: taugb1_2_6_8_12_16_lw_impl_selected = .false.
      logical :: taugb1_2_6_8_12_16_lw_entered_logged(6) = .false.
      logical :: use_native_taugb4_5_lw_impl = .false.
      logical :: taugb4_5_lw_impl_selected = .false.
      logical :: taugb4_5_lw_entered_logged(2) = .false.
      logical :: use_native_taugb3_7_9_13_15_lw_impl = .false.
      logical :: taugb3_7_9_13_15_lw_impl_selected = .false.
      logical :: taugb3_7_9_13_15_lw_entered_logged(5) = .false.
      logical :: taumol_lw_impl_logged = .false.

      contains

!----------------------------------------------------------------------------
      subroutine taumol(nlayers, pavel, wx, coldry, &
                        laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                        colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                        colbrd, fac00, fac01, fac10, fac11, &
                        rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                        rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                        rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                        selffac, selffrac, indself, forfac, forfrac, indfor, &
                        minorfrac, scaleminor, scaleminorn2, indminor, &
                        fracs, taug)
!----------------------------------------------------------------------------

      use iso_c_binding, only: c_int64_t

! *******************************************************************************
! *                                                                             *
! *                  Optical depths developed for the                           *
! *                                                                             *
! *                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
! *                                                                             *
! *                                                                             *
! *            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
! *                        131 HARTWELL AVENUE                                  *
! *                        LEXINGTON, MA 02421                                  *
! *                                                                             *
! *                                                                             *
! *                           ELI J. MLAWER                                     * 
! *                         JENNIFER DELAMERE                                   * 
! *                         STEVEN J. TAUBMAN                                   *
! *                         SHEPARD A. CLOUGH                                   *
! *                                                                             *
! *                                                                             *
! *                                                                             *
! *                                                                             *
! *                       email:  mlawer@aer.com                                *
! *                       email:  jdelamer@aer.com                              *
! *                                                                             *
! *        The authors wish to acknowledge the contributions of the             *
! *        following people:  Karen Cady-Pereira, Patrick D. Brown,             *  
! *        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
! *                                                                             *
! *******************************************************************************
! *                                                                             *
! *  Revision for g-point reduction: Michael J. Iacono, AER, Inc.               *
! *                                                                             *
! *******************************************************************************
! *     TAUMOL                                                                  *
! *                                                                             *
! *     This file contains the subroutines TAUGBn (where n goes from            *
! *     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
! *     per g-value and layer for band n.                                       *
! *                                                                             *
! *  Output:  optical depths (unitless)                                         *
! *           fractions needed to compute Planck functions at every layer       *
! *               and g-value                                                   *
! *                                                                             *
! *     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
! *     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
! *                                                                             *
! *  Input                                                                      *
! *                                                                             *
! *     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
! *     COMMON /PRECISE/  ONEMINUS                                              *
! *     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
! *     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
! *     COMMON /PROFDATA/ LAYTROP,                                              *
! *    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
! *    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),             *
! *    &                  COLO2(MXLAY)
! *     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
! *    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
! *     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
! *     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
! *                                                                             *
! *     Description:                                                            *
! *     NG(IBAND) - number of g-values in band IBAND                            *
! *     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
! *                   atmospheres that are stored for band IBAND per            *
! *                   pressure level and temperature.  Each of these            *
! *                   atmospheres has different relative amounts of the         *
! *                   key species for the band (i.e. different binary           *
! *                   species parameters).                                      *
! *     NSPB(IBAND) - same for upper atmosphere                                 *
! *     ONEMINUS - since problems are caused in some cases by interpolation     *
! *                parameters equal to or greater than 1, for these cases       *
! *                these parameters are set to this value, slightly < 1.        *
! *     PAVEL - layer pressures (mb)                                            *
! *     TAVEL - layer temperatures (degrees K)                                  *
! *     PZ - level pressures (mb)                                               *
! *     TZ - level temperatures (degrees K)                                     *
! *     LAYTROP - layer at which switch is made from one combination of         *
! *               key species to another                                        *
! *     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
! *               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
! *               respectively (molecules/cm**2)                                *
! *     FACij(LAY) - for layer LAY, these are factors that are needed to        *
! *                  compute the interpolation factors that multiply the        *
! *                  appropriate reference k-values.  A value of 0 (1) for      *
! *                  i,j indicates that the corresponding factor multiplies     *
! *                  reference k-value for the lower (higher) of the two        *
! *                  appropriate temperatures, and altitudes, respectively.     *
! *     JP - the index of the lower (in altitude) of the two appropriate        *
! *          reference pressure levels needed for interpolation                 *
! *     JT, JT1 - the indices of the lower of the two appropriate reference     *
! *               temperatures needed for interpolation (for pressure           *
! *               levels JP and JP+1, respectively)                             *
! *     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
! *               (water vapor density)/(atmospheric density at 296K and        *
! *               1013 mb)                                                      *
! *     SELFFRAC - factor needed for temperature interpolation of reference     *
! *                water vapor self-continuum data                              *
! *     INDSELF - index of the lower of the two appropriate reference           *
! *               temperatures needed for the self-continuum interpolation      *
! *     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
! *     FORFRAC - factor needed for temperature interpolation of reference      *
! *                water vapor foreign-continuum data                           *
! *     INDFOR  - index of the lower of the two appropriate reference           *
! *               temperatures needed for the foreign-continuum interpolation   *
! *                                                                             *
! *  Data input                                                                 *
! *     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
! *                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
! *        (note:  n is the band number,'MGAS' is the species name of the minor *
! *         gas)                                                                *
! *                                                                             *
! *     Description:                                                            *
! *     KA - k-values for low reference atmospheres (key-species only)          *
! *          (units: cm**2/molecule)                                            *
! *     KB - k-values for high reference atmospheres (key-species only)         *
! *          (units: cm**2/molecule)                                            *
! *     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
! *          (units: cm**2/molecule)                                            *
! *     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
! *          (units: cm**2/molecule)                                            *
! *     SELFREF - k-values for water vapor self-continuum for reference         *
! *               atmospheres (used below LAYTROP)                              *
! *               (units: cm**2/molecule)                                       *
! *     FORREF  - k-values for water vapor foreign-continuum for reference      *
! *               atmospheres (used below/above LAYTROP)                        *
! *               (units: cm**2/molecule)                                       *
! *                                                                             *
! *     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
! *     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
! *                                                                             *
!*******************************************************************************

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: nlayers         ! total number of layers
      real(kind=r8), intent(in) :: pavel(:)           ! layer pressures (mb) 
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: wx(:,:)            ! cross-section amounts (mol/cm2)
                                                      !    Dimensions: (maxxsec,nlayers)
      real(kind=r8), intent(in) :: coldry(:)          ! column amount (dry air)
                                                      !    Dimensions: (nlayers)

      integer, intent(in) :: laytrop         ! tropopause layer index
      integer, intent(in) :: jp(:)           ! 
                                                      !    Dimensions: (nlayers)
      integer, intent(in) :: jt(:)           !
                                                      !    Dimensions: (nlayers)
      integer, intent(in) :: jt1(:)          !
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: planklay(:,:)      ! 
                                                      !    Dimensions: (nlayers,nbndlw)
      real(kind=r8), intent(in) :: planklev(0:,:)     ! 
                                                      !    Dimensions: (nlayers,nbndlw)
      real(kind=r8), intent(in) :: plankbnd(:)        ! 
                                                      !    Dimensions: (nbndlw)

      real(kind=r8), intent(in) :: colh2o(:)          ! column amount (h2o)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colco2(:)          ! column amount (co2)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colo3(:)           ! column amount (o3)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: coln2o(:)          ! column amount (n2o)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colco(:)           ! column amount (co)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colch4(:)          ! column amount (ch4)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colo2(:)           ! column amount (o2)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colbrd(:)          ! column amount (broadening gases)
                                                      !    Dimensions: (nlayers)

      integer, intent(in) :: indself(:)
                                                      !    Dimensions: (nlayers)
      integer, intent(in) :: indfor(:)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: selffac(:)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: selffrac(:)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: forfac(:)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: forfrac(:)
                                                      !    Dimensions: (nlayers)

      integer, intent(in) :: indminor(:)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: minorfrac(:)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: scaleminor(:)
                                                      !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: scaleminorn2(:)
                                                      !    Dimensions: (nlayers)

      real(kind=r8), intent(in) :: &                  !
                       fac00(:), fac01(:), &          !    Dimensions: (nlayers)
                       fac10(:), fac11(:) 
      real(kind=r8), intent(in) :: &                  !
                       rat_h2oco2(:),rat_h2oco2_1(:), &
                       rat_h2oo3(:),rat_h2oo3_1(:), & !    Dimensions: (nlayers)
                       rat_h2on2o(:),rat_h2on2o_1(:), &
                       rat_h2och4(:),rat_h2och4_1(:), &
                       rat_n2oco2(:),rat_n2oco2_1(:), &
                       rat_o3co2(:),rat_o3co2_1(:)

! ----- Output -----
      real(kind=r8), intent(out) :: fracs(:,:)        ! planck fractions
                                                      !    Dimensions: (nlayers,ngptlw)
      real(kind=r8), intent(out) :: taug(:,:)         ! gaseous optical depth 
                                                      !    Dimensions: (nlayers,ngptlw)
      integer(c_int64_t), target :: jp64(nlayers), jt64(nlayers), jt164(nlayers)
      integer(c_int64_t), target :: indself64(nlayers), indfor64(nlayers), indminor64(nlayers)
      integer :: lay_idx

      hvrtau = '$Revision: 1.7 $'

      do lay_idx = 1, nlayers
         jp64(lay_idx) = int(jp(lay_idx), c_int64_t)
         jt64(lay_idx) = int(jt(lay_idx), c_int64_t)
         jt164(lay_idx) = int(jt1(lay_idx), c_int64_t)
         indself64(lay_idx) = int(indself(lay_idx), c_int64_t)
         indfor64(lay_idx) = int(indfor(lay_idx), c_int64_t)
         indminor64(lay_idx) = int(indminor(lay_idx), c_int64_t)
      enddo

! Calculate gaseous optical depth and planck fractions for each spectral band.

      call taumol_lw_log_impl()
      call taugb1
      call taugb2
      call taugb3
      call taugb4
      call taugb5
      call taugb6
      call taugb7
      call taugb8
      call taugb9
      call taugb10
      call taugb11
      call taugb12
      call taugb13
      call taugb14
      call taugb15
      call taugb16

      contains

!----------------------------------------------------------------------------
      subroutine taugb1
!----------------------------------------------------------------------------

! ------- Modifications -------
!  Written by Eli J. Mlawer, Atmospheric & Environmental Research.
!  Revised by Michael J. Iacono, Atmospheric & Environmental Research.
!
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                          (high key - h2o; high minor - n2)
!
!     note: previous versions of rrtm band 1: 
!           10-250 cm-1 (low - h2o; high - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng1
      use rrlw_kg01, only : fracrefa, fracrefb, absa, absb, &
                            ka_mn2, kb_mn2, selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      real(kind=r8) :: pp, corradj, scalen2, tauself, taufor, taun2

      interface
         subroutine rrtmg_lw_taugb1_codon(nlayers_c, laytrop_c, ng1_c, &
              nspa1_c, nspb1_c, pavel_p, colh2o_p, colbrd_p, jp_p, jt_p, &
              jt1_p, indself_p, indfor_p, indminor_p, fac00_p, fac01_p, &
              fac10_p, fac11_p, selffac_p, selffrac_p, forfac_p, forfrac_p, &
              scaleminorn2_p, minorfrac_p, fracrefa_p, fracrefb_p, absa_p, &
              absb_p, ka_mn2_p, kb_mn2_p, selfref_p, forref_p, fracs_p, &
              taug_p) bind(c, name="rrtmg_lw_taugb1_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng1_c
            integer(c_int64_t), value :: nspa1_c, nspb1_c
            type(c_ptr), value :: pavel_p, colh2o_p, colbrd_p, jp_p, jt_p, jt1_p
            type(c_ptr), value :: indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: scaleminorn2_p, minorfrac_p
            type(c_ptr), value :: fracrefa_p, fracrefb_p, absa_p, absb_p
            type(c_ptr), value :: ka_mn2_p, kb_mn2_p, selfref_p, forref_p
            type(c_ptr), value :: fracs_p, taug_p
         end subroutine rrtmg_lw_taugb1_codon
      end interface


! Minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

      call taugb1_2_6_8_12_16_lw_select_impl()
      if (.not. use_native_taugb1_2_6_8_12_16_lw_impl) then
         call taugb1_2_6_8_12_16_lw_log_entered(1)
         call rrtmg_lw_taugb1_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng1, c_int64_t), &
              int(nspa(1), c_int64_t), int(nspb(1), c_int64_t), &
              transfer(loc(pavel(1)), c_null_ptr), transfer(loc(colh2o(1)), c_null_ptr), &
              transfer(loc(colbrd(1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(scaleminorn2(1)), c_null_ptr), &
              transfer(loc(minorfrac(1)), c_null_ptr), transfer(loc(fracrefa(1)), c_null_ptr), &
              transfer(loc(fracrefb(1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(ka_mn2(1,1)), c_null_ptr), &
              transfer(loc(kb_mn2(1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(1) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(1) + 1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)
         pp = pavel(lay)
         corradj =  1.
         if (pp .lt. 250._r8) then
            corradj = 1._r8 - 0.15_r8 * (250._r8-pp) / 154.4_r8
         endif

         scalen2 = colbrd(lay) * scaleminorn2(lay)
         do ig = 1, ng1
            tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) -  forref(indf,ig))) 
            taun2 = scalen2*(ka_mn2(indm,ig) + & 
                 minorfrac(lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
            taug(lay,ig) = corradj * (colh2o(lay) * &
                (fac00(lay) * absa(ind0,ig) + &
                 fac10(lay) * absa(ind0+1,ig) + &
                 fac01(lay) * absa(ind1,ig) + &
                 fac11(lay) * absa(ind1+1,ig)) & 
                 + tauself + taufor + taun2)
             fracs(lay,ig) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(1) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(1) + 1
         indf = indfor(lay)
         indm = indminor(lay)
         pp = pavel(lay)
         corradj =  1._r8 - 0.15_r8 * (pp / 95.6_r8)

         scalen2 = colbrd(lay) * scaleminorn2(lay)
         do ig = 1, ng1
            taufor = forfac(lay) * (forref(indf,ig) + &
                 forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig))) 
            taun2 = scalen2*(kb_mn2(indm,ig) + & 
                 minorfrac(lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
            taug(lay,ig) = corradj * (colh2o(lay) * &
                (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) + &
                 fac11(lay) * absb(ind1+1,ig)) &  
                 + taufor + taun2)
            fracs(lay,ig) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb1

!----------------------------------------------------------------------------
      subroutine taugb2
!----------------------------------------------------------------------------
!
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!
!     note: previous version of rrtm band 2: 
!           250 - 500 cm-1 (low - h2o; high - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng2, ngs1
      use rrlw_kg02, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, ig
      real(kind=r8) :: pp, corradj, tauself, taufor

      interface
         subroutine rrtmg_lw_taugb2_codon(nlayers_c, laytrop_c, ng2_c, ngs1_c, &
              nspa2_c, nspb2_c, pavel_p, colh2o_p, jp_p, jt_p, jt1_p, &
              indself_p, indfor_p, fac00_p, fac01_p, fac10_p, fac11_p, &
              selffac_p, selffrac_p, forfac_p, forfrac_p, fracrefa_p, &
              fracrefb_p, absa_p, absb_p, selfref_p, forref_p, fracs_p, &
              taug_p) bind(c, name="rrtmg_lw_taugb2_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng2_c, ngs1_c
            integer(c_int64_t), value :: nspa2_c, nspb2_c
            type(c_ptr), value :: pavel_p, colh2o_p, jp_p, jt_p, jt1_p
            type(c_ptr), value :: indself_p, indfor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: fracrefa_p, fracrefb_p, absa_p, absb_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb2_codon
      end interface


! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

      call taugb1_2_6_8_12_16_lw_select_impl()
      if (.not. use_native_taugb1_2_6_8_12_16_lw_impl) then
         call taugb1_2_6_8_12_16_lw_log_entered(2)
         call rrtmg_lw_taugb2_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng2, c_int64_t), &
              int(ngs1, c_int64_t), int(nspa(2), c_int64_t), int(nspb(2), c_int64_t), &
              transfer(loc(pavel(1)), c_null_ptr), transfer(loc(colh2o(1)), c_null_ptr), &
              c_loc(jp64(1)), c_loc(jt64(1)), c_loc(jt164(1)), c_loc(indself64(1)), &
              c_loc(indfor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(fracrefa(1)), c_null_ptr), &
              transfer(loc(fracrefb(1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(2) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(2) + 1
         inds = indself(lay)
         indf = indfor(lay)
         pp = pavel(lay)
         corradj = 1._r8 - .05_r8 * (pp - 100._r8) / 900._r8
         do ig = 1, ng2
            tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs1+ig) = corradj * (colh2o(lay) * &
                (fac00(lay) * absa(ind0,ig) + &
                 fac10(lay) * absa(ind0+1,ig) + &
                 fac01(lay) * absa(ind1,ig) + &
                 fac11(lay) * absa(ind1+1,ig)) &
                 + tauself + taufor)
            fracs(lay,ngs1+ig) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(2) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(2) + 1
         indf = indfor(lay)
         do ig = 1, ng2
            taufor =  forfac(lay) * (forref(indf,ig) + &
                 forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs1+ig) = colh2o(lay) * &
                (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) + &
                 fac11(lay) * absb(ind1+1,ig)) &
                 + taufor
            fracs(lay,ngs1+ig) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb2

!----------------------------------------------------------------------------
      subroutine taugb3
!----------------------------------------------------------------------------
!
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng3, ngs2
      use rrlw_ref, only : chi_mls
      use rrlw_kg03, only : fracrefa, fracrefb, absa, absb,&
                            ka_mn2o, kb_mn2o, selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmn2o, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
                       fmn2o, fmn2omf, chi_n2o, ratn2o, adjfac, adjcoln2o
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, n2om1, n2om2, absn2o
      real(kind=r8) :: refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb3_codon(nlayers_c, laytrop_c, ng3_c, ngs2_c, &
              nspa3_c, nspb3_c, oneminus_c, colh2o_p, colco2_p, coln2o_p, &
              coldry_p, rat_h2oco2_p, rat_h2oco2_1_p, jp_p, jt_p, jt1_p, &
              indself_p, indfor_p, indminor_p, fac00_p, fac01_p, fac10_p, &
              fac11_p, selffac_p, selffrac_p, forfac_p, forfrac_p, &
              minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p, absa_p, &
              absb_p, ka_mn2o_p, kb_mn2o_p, selfref_p, forref_p, fracs_p, &
              taug_p) bind(c, name="rrtmg_lw_taugb3_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng3_c, ngs2_c
            integer(c_int64_t), value :: nspa3_c, nspb3_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: colh2o_p, colco2_p, coln2o_p, coldry_p
            type(c_ptr), value :: rat_h2oco2_p, rat_h2oco2_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p
            type(c_ptr), value :: absa_p, absb_p, ka_mn2o_p, kb_mn2o_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb3_codon
      end interface


! Minor gas mapping levels:
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

!  P = 212.725 mb
      refrat_planck_a = chi_mls(1,9)/chi_mls(2,9)

!  P = 95.58 mb
      refrat_planck_b = chi_mls(1,13)/chi_mls(2,13)

!  P = 706.270mb
      refrat_m_a = chi_mls(1,3)/chi_mls(2,3)

!  P = 95.58 mb 
      refrat_m_b = chi_mls(1,13)/chi_mls(2,13)

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  Below laytrop, the water vapor 
! self-continuum and foreign continuum is interpolated (in temperature) 
! separately.

      call taugb3_7_9_13_15_lw_select_impl()
      if (.not. use_native_taugb3_7_9_13_15_lw_impl) then
         call taugb3_7_9_13_15_lw_log_entered(3)
         call rrtmg_lw_taugb3_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng3, c_int64_t), &
              int(ngs2, c_int64_t), int(nspa(3), c_int64_t), int(nspb(3), c_int64_t), &
              real(oneminus, c_double), transfer(loc(colh2o(1)), c_null_ptr), &
              transfer(loc(colco2(1)), c_null_ptr), transfer(loc(coln2o(1)), c_null_ptr), &
              transfer(loc(coldry(1)), c_null_ptr), transfer(loc(rat_h2oco2(1)), c_null_ptr), &
              transfer(loc(rat_h2oco2_1(1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(minorfrac(1)), c_null_ptr), &
              transfer(loc(chi_mls(1,1)), c_null_ptr), transfer(loc(fracrefa(1,1)), c_null_ptr), &
              transfer(loc(fracrefb(1,1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(ka_mn2o(1,1,1)), c_null_ptr), &
              transfer(loc(kb_mn2o(1,1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)        

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mn2o = colh2o(lay) + refrat_m_a*colco2(lay)
         specparm_mn2o = colh2o(lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 8._r8*specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0_r8)
         fmn2omf = minorfrac(lay)*fmn2o
!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(lay)/coldry(lay)
         ratn2o = 1.e20_r8*chi_n2o/chi_mls(4,jp(lay)+1)
         if (ratn2o .gt. 1.5_r8) then
            adjfac = 0.5_r8+(ratn2o-0.5_r8)**0.65_r8
            adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_r8
         else
            adjcoln2o = coln2o(lay)
         endif

         speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(3) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(3) + js1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif
         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng3
            tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) +  &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs2+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(lay,ngs2+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         fac000 = (1._r8 - fs) * fac00(lay)
         fac010 = (1._r8 - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._r8 - fs1) * fac01(lay)
         fac011 = (1._r8 - fs1) * fac11(lay)
         fac101 = fs1 * fac01(lay)
         fac111 = fs1 * fac11(lay)

         speccomb_mn2o = colh2o(lay) + refrat_m_b*colco2(lay)
         specparm_mn2o = colh2o(lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 4._r8*specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0_r8)
         fmn2omf = minorfrac(lay)*fmn2o
!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(lay)/coldry(lay)
         ratn2o = 1.e20*chi_n2o/chi_mls(4,jp(lay)+1)
         if (ratn2o .gt. 1.5_r8) then
            adjfac = 0.5_r8+(ratn2o-0.5_r8)**0.65_r8
            adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_r8
         else
            adjcoln2o = coln2o(lay)
         endif

         speccomb_planck = colh2o(lay)+refrat_planck_b*colco2(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(3) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(3) + js1
         indf = indfor(lay)
         indm = indminor(lay)

         do ig = 1, ng3
            taufor = forfac(lay) * (forref(indf,ig) + &
                 forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig))) 
            n2om1 = kb_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (kb_mn2o(jmn2o+1,indm,ig)-kb_mn2o(jmn2o,indm,ig))
            n2om2 = kb_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (kb_mn2o(jmn2o+1,indm+1,ig)-kb_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)
            taug(lay,ngs2+ig) = speccomb * &
                (fac000 * absb(ind0,ig) + &
                fac100 * absb(ind0+1,ig) + &
                fac010 * absb(ind0+5,ig) + &
                fac110 * absb(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absb(ind1,ig) +  &
                fac101 * absb(ind1+1,ig) + &
                fac011 * absb(ind1+5,ig) + &
                fac111 * absb(ind1+6,ig))  &
                + taufor &
                + adjcoln2o*absn2o
            fracs(lay,ngs2+ig) = fracrefb(ig,jpl) + fpl * &
                (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
         enddo
      enddo

      end subroutine taugb3

!----------------------------------------------------------------------------
      subroutine taugb4
!----------------------------------------------------------------------------
!
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng4, ngs3
      use rrlw_ref, only : chi_mls
      use rrlw_kg04, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor
      real(kind=r8) :: refrat_planck_a, refrat_planck_b
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb4_codon(nlayers_c, laytrop_c, ng4_c, ngs3_c, &
              nspa4_c, nspb4_c, oneminus_c, colh2o_p, colco2_p, colo3_p, &
              rat_h2oco2_p, rat_h2oco2_1_p, rat_o3co2_p, rat_o3co2_1_p, &
              jp_p, jt_p, jt1_p, indself_p, indfor_p, fac00_p, fac01_p, &
              fac10_p, fac11_p, selffac_p, selffrac_p, forfac_p, forfrac_p, &
              chi_mls_p, fracrefa_p, fracrefb_p, absa_p, absb_p, selfref_p, &
              forref_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb4_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng4_c, ngs3_c
            integer(c_int64_t), value :: nspa4_c, nspb4_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: colh2o_p, colco2_p, colo3_p
            type(c_ptr), value :: rat_h2oco2_p, rat_h2oco2_1_p, rat_o3co2_p, rat_o3co2_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: chi_mls_p, fracrefa_p, fracrefb_p, absa_p, absb_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb4_codon
      end interface


! P =   142.5940 mb
      refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)

! P = 95.58350 mb
      refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  Below laytrop, the water 
! vapor self-continuum and foreign continuum is interpolated (in temperature) 
! separately.

      call taugb4_5_lw_select_impl()
      if (.not. use_native_taugb4_5_lw_impl) then
         call taugb4_5_lw_log_entered(4)
         call rrtmg_lw_taugb4_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng4, c_int64_t), &
              int(ngs3, c_int64_t), int(nspa(4), c_int64_t), int(nspb(4), c_int64_t), &
              real(oneminus, c_double), transfer(loc(colh2o(1)), c_null_ptr), &
              transfer(loc(colco2(1)), c_null_ptr), transfer(loc(colo3(1)), c_null_ptr), &
              transfer(loc(rat_h2oco2(1)), c_null_ptr), transfer(loc(rat_h2oco2_1(1)), c_null_ptr), &
              transfer(loc(rat_o3co2(1)), c_null_ptr), transfer(loc(rat_o3co2_1(1)), c_null_ptr), &
              c_loc(jp64(1)), c_loc(jt64(1)), c_loc(jt164(1)), c_loc(indself64(1)), &
              c_loc(indfor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(chi_mls(1,1)), c_null_ptr), &
              transfer(loc(fracrefa(1,1)), c_null_ptr), transfer(loc(fracrefb(1,1)), c_null_ptr), &
              transfer(loc(absa(1,1)), c_null_ptr), transfer(loc(absb(1,1)), c_null_ptr), &
              transfer(loc(selfref(1,1)), c_null_ptr), transfer(loc(forref(1,1)), c_null_ptr), &
              transfer(loc(fracs(1,1)), c_null_ptr), transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(4) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(4) + js1
         inds = indself(lay)
         indf = indfor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng4
            tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) +  &
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs3+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracs(lay,ngs3+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         speccomb = colo3(lay) + rat_o3co2(lay)*colco2(lay)
         specparm = colo3(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colo3(lay) + rat_o3co2_1(lay)*colco2(lay)
         specparm1 = colo3(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         fac000 = (1._r8 - fs) * fac00(lay)
         fac010 = (1._r8 - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._r8 - fs1) * fac01(lay)
         fac011 = (1._r8 - fs1) * fac11(lay)
         fac101 = fs1 * fac01(lay)
         fac111 = fs1 * fac11(lay)

         speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay)
         specparm_planck = colo3(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(4) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(4) + js1

         do ig = 1, ng4
            taug(lay,ngs3+ig) =  speccomb * &
                (fac000 * absb(ind0,ig) + &
                fac100 * absb(ind0+1,ig) + &
                fac010 * absb(ind0+5,ig) + &
                fac110 * absb(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absb(ind1,ig) +  &
                fac101 * absb(ind1+1,ig) + &
                fac011 * absb(ind1+5,ig) + &
                fac111 * absb(ind1+6,ig))
            fracs(lay,ngs3+ig) = fracrefb(ig,jpl) + fpl * &
                (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for co2.  Revised to apply weighting for g-point reduction in this band.

         taug(lay,ngs3+8)=taug(lay,ngs3+8)*0.92
         taug(lay,ngs3+9)=taug(lay,ngs3+9)*0.88
         taug(lay,ngs3+10)=taug(lay,ngs3+10)*1.07
         taug(lay,ngs3+11)=taug(lay,ngs3+11)*1.1
         taug(lay,ngs3+12)=taug(lay,ngs3+12)*0.99
         taug(lay,ngs3+13)=taug(lay,ngs3+13)*0.88
         taug(lay,ngs3+14)=taug(lay,ngs3+14)*0.943

      enddo

      end subroutine taugb4

!----------------------------------------------------------------------------
      subroutine taugb5
!----------------------------------------------------------------------------
!
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng5, ngs4
      use rrlw_ref, only : chi_mls
      use rrlw_kg05, only : fracrefa, fracrefb, absa, absb, &
                            ka_mo3, selfref, forref, ccl4

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmo3, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mo3, specparm_mo3, specmult_mo3, fmo3
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, o3m1, o3m2, abso3
      real(kind=r8) :: refrat_planck_a, refrat_planck_b, refrat_m_a
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb5_codon(nlayers_c, laytrop_c, ng5_c, ngs4_c, &
              nspa5_c, nspb5_c, maxxsec_c, oneminus_c, colh2o_p, colco2_p, &
              colo3_p, wx_p, rat_h2oco2_p, rat_h2oco2_1_p, rat_o3co2_p, &
              rat_o3co2_1_p, jp_p, jt_p, jt1_p, indself_p, indfor_p, &
              indminor_p, fac00_p, fac01_p, fac10_p, fac11_p, selffac_p, &
              selffrac_p, forfac_p, forfrac_p, minorfrac_p, chi_mls_p, &
              fracrefa_p, fracrefb_p, absa_p, absb_p, ka_mo3_p, selfref_p, &
              forref_p, ccl4_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb5_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng5_c, ngs4_c
            integer(c_int64_t), value :: nspa5_c, nspb5_c, maxxsec_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: colh2o_p, colco2_p, colo3_p, wx_p
            type(c_ptr), value :: rat_h2oco2_p, rat_h2oco2_1_p, rat_o3co2_p, rat_o3co2_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p, minorfrac_p
            type(c_ptr), value :: chi_mls_p, fracrefa_p, fracrefb_p, absa_p, absb_p
            type(c_ptr), value :: ka_mo3_p, selfref_p, forref_p, ccl4_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb5_codon
      end interface


! Minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

! P = 473.420 mb
      refrat_planck_a = chi_mls(1,5)/chi_mls(2,5)

! P = 0.2369 mb
      refrat_planck_b = chi_mls(3,43)/chi_mls(2,43)

! P = 317.3480
      refrat_m_a = chi_mls(1,7)/chi_mls(2,7)

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  Below laytrop, the 
! water vapor self-continuum and foreign continuum is 
! interpolated (in temperature) separately.

      call taugb4_5_lw_select_impl()
      if (.not. use_native_taugb4_5_lw_impl) then
         call taugb4_5_lw_log_entered(5)
         call rrtmg_lw_taugb5_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng5, c_int64_t), &
              int(ngs4, c_int64_t), int(nspa(5), c_int64_t), int(nspb(5), c_int64_t), &
              int(maxxsec, c_int64_t), real(oneminus, c_double), &
              transfer(loc(colh2o(1)), c_null_ptr), transfer(loc(colco2(1)), c_null_ptr), &
              transfer(loc(colo3(1)), c_null_ptr), transfer(loc(wx(1,1)), c_null_ptr), &
              transfer(loc(rat_h2oco2(1)), c_null_ptr), transfer(loc(rat_h2oco2_1(1)), c_null_ptr), &
              transfer(loc(rat_o3co2(1)), c_null_ptr), transfer(loc(rat_o3co2_1(1)), c_null_ptr), &
              c_loc(jp64(1)), c_loc(jt64(1)), c_loc(jt164(1)), c_loc(indself64(1)), &
              c_loc(indfor64(1)), c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(minorfrac(1)), c_null_ptr), &
              transfer(loc(chi_mls(1,1)), c_null_ptr), transfer(loc(fracrefa(1,1)), c_null_ptr), &
              transfer(loc(fracrefb(1,1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(ka_mo3(1,1,1)), c_null_ptr), &
              transfer(loc(selfref(1,1)), c_null_ptr), transfer(loc(forref(1,1)), c_null_ptr), &
              transfer(loc(ccl4(1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mo3 = colh2o(lay) + refrat_m_a*colco2(lay)
         specparm_mo3 = colh2o(lay)/speccomb_mo3
         if (specparm_mo3 .ge. oneminus) specparm_mo3 = oneminus
         specmult_mo3 = 8._r8*specparm_mo3
         jmo3 = 1 + int(specmult_mo3)
         fmo3 = mod(specmult_mo3,1.0_r8)

         speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(5) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(5) + js1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng5
            tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            o3m1 = ka_mo3(jmo3,indm,ig) + fmo3 * &
                 (ka_mo3(jmo3+1,indm,ig)-ka_mo3(jmo3,indm,ig))
            o3m2 = ka_mo3(jmo3,indm+1,ig) + fmo3 * &
                 (ka_mo3(jmo3+1,indm+1,ig)-ka_mo3(jmo3,indm+1,ig))
            abso3 = o3m1 + minorfrac(lay)*(o3m2-o3m1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * & 
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs4+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + abso3*colo3(lay) &
                 + wx(1,lay) * ccl4(ig)
            fracs(lay,ngs4+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

         speccomb = colo3(lay) + rat_o3co2(lay)*colco2(lay)
         specparm = colo3(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colo3(lay) + rat_o3co2_1(lay)*colco2(lay)
         specparm1 = colo3(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         fac000 = (1._r8 - fs) * fac00(lay)
         fac010 = (1._r8 - fs) * fac10(lay)
         fac100 = fs * fac00(lay)
         fac110 = fs * fac10(lay)
         fac001 = (1._r8 - fs1) * fac01(lay)
         fac011 = (1._r8 - fs1) * fac11(lay)
         fac101 = fs1 * fac01(lay)
         fac111 = fs1 * fac11(lay)

         speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay)
         specparm_planck = colo3(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(5) + js
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(5) + js1
         
         do ig = 1, ng5
            taug(lay,ngs4+ig) = speccomb * &
                (fac000 * absb(ind0,ig) + &
                fac100 * absb(ind0+1,ig) + &
                fac010 * absb(ind0+5,ig) + &
                fac110 * absb(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absb(ind1,ig) + &
                fac101 * absb(ind1+1,ig) + &
                fac011 * absb(ind1+5,ig) + &
                fac111 * absb(ind1+6,ig))  &
                + wx(1,lay) * ccl4(ig)
            fracs(lay,ngs4+ig) = fracrefb(ig,jpl) + fpl * &
                (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
         enddo
      enddo

      end subroutine taugb5

!----------------------------------------------------------------------------
      subroutine taugb6
!----------------------------------------------------------------------------
!
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng6, ngs5
      use rrlw_ref, only : chi_mls
      use rrlw_kg06, only : fracrefa, absa, ka_mco2, &
                            selfref, forref, cfc11adj, cfc12

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      real(kind=r8) :: chi_co2, ratco2, adjfac, adjcolco2
      real(kind=r8) :: tauself, taufor, absco2

      interface
         subroutine rrtmg_lw_taugb6_codon(nlayers_c, laytrop_c, ng6_c, ngs5_c, &
              nspa6_c, maxxsec_c, colh2o_p, colco2_p, coldry_p, wx_p, jp_p, jt_p, &
              jt1_p, indself_p, indfor_p, indminor_p, fac00_p, fac01_p, fac10_p, &
              fac11_p, selffac_p, selffrac_p, forfac_p, forfrac_p, minorfrac_p, &
              chi_mls_p, fracrefa_p, absa_p, ka_mco2_p, selfref_p, forref_p, &
              cfc11adj_p, cfc12_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb6_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng6_c, ngs5_c
            integer(c_int64_t), value :: nspa6_c, maxxsec_c
            type(c_ptr), value :: colh2o_p, colco2_p, coldry_p, wx_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: minorfrac_p, chi_mls_p, fracrefa_p, absa_p
            type(c_ptr), value :: ka_mco2_p, selfref_p, forref_p
            type(c_ptr), value :: cfc11adj_p, cfc12_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb6_codon
      end interface


! Minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12

! Compute the optical depth by interpolating in ln(pressure) and
! temperature. The water vapor self-continuum and foreign continuum
! is interpolated (in temperature) separately.  

      call taugb1_2_6_8_12_16_lw_select_impl()
      if (.not. use_native_taugb1_2_6_8_12_16_lw_impl) then
         call taugb1_2_6_8_12_16_lw_log_entered(6)
         call rrtmg_lw_taugb6_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng6, c_int64_t), &
              int(ngs5, c_int64_t), int(nspa(6), c_int64_t), int(maxxsec, c_int64_t), &
              transfer(loc(colh2o(1)), c_null_ptr), transfer(loc(colco2(1)), c_null_ptr), &
              transfer(loc(coldry(1)), c_null_ptr), transfer(loc(wx(1,1)), c_null_ptr), &
              c_loc(jp64(1)), c_loc(jt64(1)), c_loc(jt164(1)), c_loc(indself64(1)), &
              c_loc(indfor64(1)), c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(minorfrac(1)), c_null_ptr), &
              transfer(loc(chi_mls(1,1)), c_null_ptr), transfer(loc(fracrefa(1)), c_null_ptr), &
              transfer(loc(absa(1,1)), c_null_ptr), transfer(loc(ka_mco2(1,1)), c_null_ptr), &
              transfer(loc(selfref(1,1)), c_null_ptr), transfer(loc(forref(1,1)), c_null_ptr), &
              transfer(loc(cfc11adj(1)), c_null_ptr), transfer(loc(cfc12(1)), c_null_ptr), &
              transfer(loc(fracs(1,1)), c_null_ptr), transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

! In atmospheres where the amount of CO2 is too great to be considered
! a minor species, adjust the column amount of CO2 by an empirical factor 
! to obtain the proper contribution.
         chi_co2 = colco2(lay)/(coldry(lay))
         ratco2 = 1.e20_r8*chi_co2/chi_mls(2,jp(lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.77_r8
            adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_r8
         else
            adjcolco2 = colco2(lay)
         endif

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(6) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(6) + 1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)

         do ig = 1, ng6
            tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            taug(lay,ngs5+ig) = colh2o(lay) * &
                (fac00(lay) * absa(ind0,ig) + &
                 fac10(lay) * absa(ind0+1,ig) + &
                 fac01(lay) * absa(ind1,ig) +  &
                 fac11(lay) * absa(ind1+1,ig))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx(2,lay) * cfc11adj(ig) &
                 + wx(3,lay) * cfc12(ig)
            fracs(lay,ngs5+ig) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
! Nothing important goes on above laytrop in this band.
      do lay = laytrop+1, nlayers

         do ig = 1, ng6
            taug(lay,ngs5+ig) = 0.0_r8 &
                 + wx(2,lay) * cfc11adj(ig) &
                 + wx(3,lay) * cfc12(ig)
            fracs(lay,ngs5+ig) = fracrefa(ig)
         enddo
      enddo

      end subroutine taugb6

!----------------------------------------------------------------------------
      subroutine taugb7
!----------------------------------------------------------------------------
!
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng7, ngs6
      use rrlw_ref, only : chi_mls
      use rrlw_kg07, only : fracrefa, fracrefb, absa, absb, &
                            ka_mco2, kb_mco2, selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmco2, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, co2m1, co2m2, absco2
      real(kind=r8) :: chi_co2, ratco2, adjfac, adjcolco2
      real(kind=r8) :: refrat_planck_a, refrat_m_a
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb7_codon(nlayers_c, laytrop_c, ng7_c, ngs6_c, &
              nspa7_c, nspb7_c, oneminus_c, colh2o_p, colco2_p, colo3_p, &
              coldry_p, rat_h2oo3_p, rat_h2oo3_1_p, jp_p, jt_p, jt1_p, &
              indself_p, indfor_p, indminor_p, fac00_p, fac01_p, fac10_p, &
              fac11_p, selffac_p, selffrac_p, forfac_p, forfrac_p, &
              minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p, absa_p, &
              absb_p, ka_mco2_p, kb_mco2_p, selfref_p, forref_p, fracs_p, &
              taug_p) bind(c, name="rrtmg_lw_taugb7_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng7_c, ngs6_c
            integer(c_int64_t), value :: nspa7_c, nspb7_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: colh2o_p, colco2_p, colo3_p, coldry_p
            type(c_ptr), value :: rat_h2oo3_p, rat_h2oo3_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p
            type(c_ptr), value :: absa_p, absb_p, ka_mco2_p, kb_mco2_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb7_codon
      end interface


! Minor gas mapping level :
!     lower - co2, p = 706.2620 mbar, t= 278.94 k
!     upper - co2, p = 12.9350 mbar, t = 234.01 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.

! P = 706.2620 mb
      refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)

! P = 706.2720 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(3,3)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately. 

      call taugb3_7_9_13_15_lw_select_impl()
      if (.not. use_native_taugb3_7_9_13_15_lw_impl) then
         call taugb3_7_9_13_15_lw_log_entered(7)
         call rrtmg_lw_taugb7_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng7, c_int64_t), &
              int(ngs6, c_int64_t), int(nspa(7), c_int64_t), int(nspb(7), c_int64_t), &
              real(oneminus, c_double), transfer(loc(colh2o(1)), c_null_ptr), &
              transfer(loc(colco2(1)), c_null_ptr), transfer(loc(colo3(1)), c_null_ptr), &
              transfer(loc(coldry(1)), c_null_ptr), transfer(loc(rat_h2oo3(1)), c_null_ptr), &
              transfer(loc(rat_h2oo3_1(1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(minorfrac(1)), c_null_ptr), &
              transfer(loc(chi_mls(1,1)), c_null_ptr), transfer(loc(fracrefa(1,1)), c_null_ptr), &
              transfer(loc(fracrefb(1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(ka_mco2(1,1,1)), c_null_ptr), &
              transfer(loc(kb_mco2(1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oo3(lay)*colo3(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(lay) + rat_h2oo3_1(lay)*colo3(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mco2 = colh2o(lay) + refrat_m_a*colo3(lay)
         specparm_mco2 = colh2o(lay)/speccomb_mco2
         if (specparm_mco2 .ge. oneminus) specparm_mco2 = oneminus
         specmult_mco2 = 8._r8*specparm_mco2

         jmco2 = 1 + int(specmult_mco2)
         fmco2 = mod(specmult_mco2,1.0_r8)

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay)/(coldry(lay))
         ratco2 = 1.e20*chi_co2/chi_mls(2,jp(lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 3.0_r8+(ratco2-3.0_r8)**0.79_r8
            adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_r8
         else
            adjcolco2 = colco2(lay)
         endif

         speccomb_planck = colh2o(lay)+refrat_planck_a*colo3(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(7) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(7) + js1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif
         if (specparm .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng7
            tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) +  &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs6+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcolco2*absco2
            fracs(lay,ngs6+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay)/(coldry(lay))
         ratco2 = 1.e20*chi_co2/chi_mls(2,jp(lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.79_r8
            adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_r8
         else
            adjcolco2 = colco2(lay)
         endif

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(7) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(7) + 1
         indm = indminor(lay)

         do ig = 1, ng7
            absco2 = kb_mco2(indm,ig) + minorfrac(lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig))
            taug(lay,ngs6+ig) = colo3(lay) * &
                 (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) + &
                 fac11(lay) * absb(ind1+1,ig)) &
                 + adjcolco2 * absco2
            fracs(lay,ngs6+ig) = fracrefb(ig)
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for o3.  Revised to apply weighting for g-point reduction in this band.

         taug(lay,ngs6+6)=taug(lay,ngs6+6)*0.92_r8
         taug(lay,ngs6+7)=taug(lay,ngs6+7)*0.88_r8
         taug(lay,ngs6+8)=taug(lay,ngs6+8)*1.07_r8
         taug(lay,ngs6+9)=taug(lay,ngs6+9)*1.1_r8
         taug(lay,ngs6+10)=taug(lay,ngs6+10)*0.99_r8
         taug(lay,ngs6+11)=taug(lay,ngs6+11)*0.855_r8

      enddo

      end subroutine taugb7

!----------------------------------------------------------------------------
      subroutine taugb8
!----------------------------------------------------------------------------
!
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng8, ngs7
      use rrlw_ref, only : chi_mls
      use rrlw_kg08, only : fracrefa, fracrefb, absa, absb, &
                            ka_mco2, ka_mn2o, ka_mo3, kb_mco2, kb_mn2o, &
                            selfref, forref, cfc12, cfc22adj

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      real(kind=r8) :: tauself, taufor, absco2, abso3, absn2o
      real(kind=r8) :: chi_co2, ratco2, adjfac, adjcolco2

      interface
         subroutine rrtmg_lw_taugb8_codon(nlayers_c, laytrop_c, ng8_c, ngs7_c, &
              nspa8_c, nspb8_c, maxxsec_c, colh2o_p, colco2_p, colo3_p, &
              coln2o_p, coldry_p, wx_p, jp_p, jt_p, jt1_p, indself_p, &
              indfor_p, indminor_p, fac00_p, fac01_p, fac10_p, fac11_p, &
              selffac_p, selffrac_p, forfac_p, forfrac_p, minorfrac_p, &
              chi_mls_p, fracrefa_p, fracrefb_p, absa_p, absb_p, ka_mco2_p, &
              ka_mn2o_p, ka_mo3_p, kb_mco2_p, kb_mn2o_p, selfref_p, forref_p, &
              cfc12_p, cfc22adj_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb8_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng8_c, ngs7_c
            integer(c_int64_t), value :: nspa8_c, nspb8_c, maxxsec_c
            type(c_ptr), value :: colh2o_p, colco2_p, colo3_p, coln2o_p, coldry_p, wx_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p
            type(c_ptr), value :: absa_p, absb_p, ka_mco2_p, ka_mn2o_p, ka_mo3_p
            type(c_ptr), value :: kb_mco2_p, kb_mn2o_p, selfref_p, forref_p
            type(c_ptr), value :: cfc12_p, cfc22adj_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb8_codon
      end interface


! Minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  Below laytrop, the water vapor 
! self-continuum and foreign continuum is interpolated (in temperature) 
! separately.

      call taugb1_2_6_8_12_16_lw_select_impl()
      if (.not. use_native_taugb1_2_6_8_12_16_lw_impl) then
         call taugb1_2_6_8_12_16_lw_log_entered(8)
         call rrtmg_lw_taugb8_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng8, c_int64_t), &
              int(ngs7, c_int64_t), int(nspa(8), c_int64_t), int(nspb(8), c_int64_t), &
              int(maxxsec, c_int64_t), transfer(loc(colh2o(1)), c_null_ptr), &
              transfer(loc(colco2(1)), c_null_ptr), transfer(loc(colo3(1)), c_null_ptr), &
              transfer(loc(coln2o(1)), c_null_ptr), transfer(loc(coldry(1)), c_null_ptr), &
              transfer(loc(wx(1,1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(minorfrac(1)), c_null_ptr), &
              transfer(loc(chi_mls(1,1)), c_null_ptr), transfer(loc(fracrefa(1)), c_null_ptr), &
              transfer(loc(fracrefb(1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(ka_mco2(1,1)), c_null_ptr), &
              transfer(loc(ka_mn2o(1,1)), c_null_ptr), transfer(loc(ka_mo3(1,1)), c_null_ptr), &
              transfer(loc(kb_mco2(1,1)), c_null_ptr), transfer(loc(kb_mn2o(1,1)), c_null_ptr), &
              transfer(loc(selfref(1,1)), c_null_ptr), transfer(loc(forref(1,1)), c_null_ptr), &
              transfer(loc(cfc12(1)), c_null_ptr), transfer(loc(cfc22adj(1)), c_null_ptr), &
              transfer(loc(fracs(1,1)), c_null_ptr), transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay)/(coldry(lay))
         ratco2 = 1.e20_r8*chi_co2/chi_mls(2,jp(lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.65_r8
            adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_r8
         else
            adjcolco2 = colco2(lay)
         endif

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(8) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(8) + 1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)

         do ig = 1, ng8
            tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            abso3 =  (ka_mo3(indm,ig) + minorfrac(lay) * &
                 (ka_mo3(indm+1,ig) - ka_mo3(indm,ig)))
            absn2o =  (ka_mn2o(indm,ig) + minorfrac(lay) * &
                 (ka_mn2o(indm+1,ig) - ka_mn2o(indm,ig)))
            taug(lay,ngs7+ig) = colh2o(lay) * &
                 (fac00(lay) * absa(ind0,ig) + &
                 fac10(lay) * absa(ind0+1,ig) + &
                 fac01(lay) * absa(ind1,ig) +  &
                 fac11(lay) * absa(ind1+1,ig)) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colo3(lay) * abso3 &
                 + coln2o(lay) * absn2o &
                 + wx(3,lay) * cfc12(ig) &
                 + wx(4,lay) * cfc22adj(ig)
            fracs(lay,ngs7+ig) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay)/coldry(lay)
         ratco2 = 1.e20_r8*chi_co2/chi_mls(2,jp(lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.65_r8
            adjcolco2 = adjfac*chi_mls(2,jp(lay)+1) * coldry(lay)*1.e-20_r8
         else
            adjcolco2 = colco2(lay)
         endif

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(8) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(8) + 1
         indm = indminor(lay)

         do ig = 1, ng8
            absco2 =  (kb_mco2(indm,ig) + minorfrac(lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
            absn2o =  (kb_mn2o(indm,ig) + minorfrac(lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
            taug(lay,ngs7+ig) = colo3(lay) * &
                 (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) + &
                 fac11(lay) * absb(ind1+1,ig)) &
                 + adjcolco2*absco2 &
                 + coln2o(lay)*absn2o & 
                 + wx(3,lay) * cfc12(ig) &
                 + wx(4,lay) * cfc22adj(ig)
            fracs(lay,ngs7+ig) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb8

!----------------------------------------------------------------------------
      subroutine taugb9
!----------------------------------------------------------------------------
!
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng9, ngs8
      use rrlw_ref, only : chi_mls
      use rrlw_kg09, only : fracrefa, fracrefb, absa, absb, &
                            ka_mn2o, kb_mn2o, selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmn2o, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, n2om1, n2om2, absn2o
      real(kind=r8) :: chi_n2o, ratn2o, adjfac, adjcoln2o
      real(kind=r8) :: refrat_planck_a, refrat_m_a
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb9_codon(nlayers_c, laytrop_c, ng9_c, ngs8_c, &
              nspa9_c, nspb9_c, oneminus_c, colh2o_p, colch4_p, coln2o_p, &
              coldry_p, rat_h2och4_p, rat_h2och4_1_p, jp_p, jt_p, jt1_p, &
              indself_p, indfor_p, indminor_p, fac00_p, fac01_p, fac10_p, &
              fac11_p, selffac_p, selffrac_p, forfac_p, forfrac_p, &
              minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p, absa_p, &
              absb_p, ka_mn2o_p, kb_mn2o_p, selfref_p, forref_p, fracs_p, &
              taug_p) bind(c, name="rrtmg_lw_taugb9_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng9_c, ngs8_c
            integer(c_int64_t), value :: nspa9_c, nspb9_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: colh2o_p, colch4_p, coln2o_p, coldry_p
            type(c_ptr), value :: rat_h2och4_p, rat_h2och4_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p
            type(c_ptr), value :: absa_p, absb_p, ka_mn2o_p, kb_mn2o_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb9_codon
      end interface


! Minor gas mapping level :
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

! P = 212 mb
      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)

! P = 706.272 mb 
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

      call taugb3_7_9_13_15_lw_select_impl()
      if (.not. use_native_taugb3_7_9_13_15_lw_impl) then
         call taugb3_7_9_13_15_lw_log_entered(9)
         call rrtmg_lw_taugb9_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng9, c_int64_t), &
              int(ngs8, c_int64_t), int(nspa(9), c_int64_t), int(nspb(9), c_int64_t), &
              real(oneminus, c_double), transfer(loc(colh2o(1)), c_null_ptr), &
              transfer(loc(colch4(1)), c_null_ptr), transfer(loc(coln2o(1)), c_null_ptr), &
              transfer(loc(coldry(1)), c_null_ptr), transfer(loc(rat_h2och4(1)), c_null_ptr), &
              transfer(loc(rat_h2och4_1(1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(minorfrac(1)), c_null_ptr), &
              transfer(loc(chi_mls(1,1)), c_null_ptr), transfer(loc(fracrefa(1,1)), c_null_ptr), &
              transfer(loc(fracrefb(1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(ka_mn2o(1,1,1)), c_null_ptr), &
              transfer(loc(kb_mn2o(1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2och4(lay)*colch4(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(lay) + rat_h2och4_1(lay)*colch4(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mn2o = colh2o(lay) + refrat_m_a*colch4(lay)
         specparm_mn2o = colh2o(lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 8._r8*specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0_r8)

!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(lay)/(coldry(lay))
         ratn2o = 1.e20_r8*chi_n2o/chi_mls(4,jp(lay)+1)
         if (ratn2o .gt. 1.5_r8) then
            adjfac = 0.5_r8+(ratn2o-0.5_r8)**0.65_r8
            adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_r8
         else
            adjcoln2o = coln2o(lay)
         endif

         speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(9) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(9) + js1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng9
            tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + & 
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs8+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(lay,ngs8+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers

!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor 
!  to obtain the proper contribution.
         chi_n2o = coln2o(lay)/(coldry(lay))
         ratn2o = 1.e20_r8*chi_n2o/chi_mls(4,jp(lay)+1)
         if (ratn2o .gt. 1.5_r8) then
            adjfac = 0.5_r8+(ratn2o-0.5_r8)**0.65_r8
            adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_r8
         else
            adjcoln2o = coln2o(lay)
         endif

         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(9) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(9) + 1
         indm = indminor(lay)

         do ig = 1, ng9
            absn2o = kb_mn2o(indm,ig) + minorfrac(lay) * &
                (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
            taug(lay,ngs8+ig) = colch4(lay) * &
                 (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) +  &
                 fac11(lay) * absb(ind1+1,ig)) &
                 + adjcoln2o*absn2o
            fracs(lay,ngs8+ig) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb9

!----------------------------------------------------------------------------
      subroutine taugb10
!----------------------------------------------------------------------------
!
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng10, ngs9
      use rrlw_kg10, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, ig
      real(kind=r8) :: tauself, taufor

      interface
         subroutine rrtmg_lw_taugb10_codon(nlayers_c, laytrop_c, ng10_c, ngs9_c, &
              nspa10_c, nspb10_c, colh2o_p, jp_p, jt_p, jt1_p, indself_p, &
              indfor_p, fac00_p, fac01_p, fac10_p, fac11_p, selffac_p, &
              selffrac_p, forfac_p, forfrac_p, fracrefa_p, fracrefb_p, absa_p, &
              absb_p, selfref_p, forref_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb10_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng10_c, ngs9_c
            integer(c_int64_t), value :: nspa10_c, nspb10_c
            type(c_ptr), value :: colh2o_p, jp_p, jt_p, jt1_p, indself_p, indfor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p, selffac_p, selffrac_p
            type(c_ptr), value :: forfac_p, forfrac_p, fracrefa_p, fracrefb_p, absa_p, absb_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb10_codon
      end interface

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

      call taugb10_11_14_lw_select_impl()
      if (.not. use_native_taugb10_11_14_lw_impl) then
         call taugb10_11_14_lw_log_entered(10)
         call rrtmg_lw_taugb10_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng10, c_int64_t), &
              int(ngs9, c_int64_t), int(nspa(10), c_int64_t), int(nspb(10), c_int64_t), &
              transfer(loc(colh2o(1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              transfer(loc(fac00(1)), c_null_ptr), transfer(loc(fac01(1)), c_null_ptr), &
              transfer(loc(fac10(1)), c_null_ptr), transfer(loc(fac11(1)), c_null_ptr), &
              transfer(loc(selffac(1)), c_null_ptr), transfer(loc(selffrac(1)), c_null_ptr), &
              transfer(loc(forfac(1)), c_null_ptr), transfer(loc(forfrac(1)), c_null_ptr), &
              transfer(loc(fracrefa(1)), c_null_ptr), transfer(loc(fracrefb(1)), c_null_ptr), &
              transfer(loc(absa(1,1)), c_null_ptr), transfer(loc(absb(1,1)), c_null_ptr), &
              transfer(loc(selfref(1,1)), c_null_ptr), transfer(loc(forref(1,1)), c_null_ptr), &
              transfer(loc(fracs(1,1)), c_null_ptr), transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(10) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(10) + 1
         inds = indself(lay)
         indf = indfor(lay)

         do ig = 1, ng10
            tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs9+ig) = colh2o(lay) * &
                 (fac00(lay) * absa(ind0,ig) + &
                 fac10(lay) * absa(ind0+1,ig) + &
                 fac01(lay) * absa(ind1,ig) + &
                 fac11(lay) * absa(ind1+1,ig))  &
                 + tauself + taufor
            fracs(lay,ngs9+ig) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(10) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(10) + 1
         indf = indfor(lay)

         do ig = 1, ng10
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs9+ig) = colh2o(lay) * &
                 (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) +  &
                 fac11(lay) * absb(ind1+1,ig)) &
                 + taufor
            fracs(lay,ngs9+ig) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb10

!----------------------------------------------------------------------------
      subroutine taugb11
!----------------------------------------------------------------------------
!
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng11, ngs10
      use rrlw_kg11, only : fracrefa, fracrefb, absa, absb, &
                            ka_mo2, kb_mo2, selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      real(kind=r8) :: scaleo2, tauself, taufor, tauo2

      interface
         subroutine rrtmg_lw_taugb11_codon(nlayers_c, laytrop_c, ng11_c, ngs10_c, &
              nspa11_c, nspb11_c, colh2o_p, colo2_p, jp_p, jt_p, jt1_p, &
              indself_p, indfor_p, indminor_p, fac00_p, fac01_p, fac10_p, &
              fac11_p, selffac_p, selffrac_p, forfac_p, forfrac_p, scaleminor_p, &
              minorfrac_p, fracrefa_p, fracrefb_p, absa_p, absb_p, ka_mo2_p, &
              kb_mo2_p, selfref_p, forref_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb11_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng11_c, ngs10_c
            integer(c_int64_t), value :: nspa11_c, nspb11_c
            type(c_ptr), value :: colh2o_p, colo2_p, jp_p, jt_p, jt1_p
            type(c_ptr), value :: indself_p, indfor_p, indminor_p, fac00_p, fac01_p
            type(c_ptr), value :: fac10_p, fac11_p, selffac_p, selffrac_p, forfac_p
            type(c_ptr), value :: forfrac_p, scaleminor_p, minorfrac_p, fracrefa_p, fracrefb_p
            type(c_ptr), value :: absa_p, absb_p, ka_mo2_p, kb_mo2_p, selfref_p, forref_p
            type(c_ptr), value :: fracs_p, taug_p
         end subroutine rrtmg_lw_taugb11_codon
      end interface

! Minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum and
! foreign continuum is interpolated (in temperature) separately.

      call taugb10_11_14_lw_select_impl()
      if (.not. use_native_taugb10_11_14_lw_impl) then
         call taugb10_11_14_lw_log_entered(11)
         call rrtmg_lw_taugb11_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng11, c_int64_t), &
              int(ngs10, c_int64_t), int(nspa(11), c_int64_t), int(nspb(11), c_int64_t), &
              transfer(loc(colh2o(1)), c_null_ptr), transfer(loc(colo2(1)), c_null_ptr), &
              c_loc(jp64(1)), c_loc(jt64(1)), c_loc(jt164(1)), c_loc(indself64(1)), &
              c_loc(indfor64(1)), c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(scaleminor(1)), c_null_ptr), &
              transfer(loc(minorfrac(1)), c_null_ptr), transfer(loc(fracrefa(1)), c_null_ptr), &
              transfer(loc(fracrefb(1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(ka_mo2(1,1)), c_null_ptr), &
              transfer(loc(kb_mo2(1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(11) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(11) + 1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)
         scaleo2 = colo2(lay)*scaleminor(lay)
         do ig = 1, ng11
            tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (ka_mo2(indm,ig) + minorfrac(lay) * &
                 (ka_mo2(indm+1,ig) - ka_mo2(indm,ig)))
            taug(lay,ngs10+ig) = colh2o(lay) * &
                 (fac00(lay) * absa(ind0,ig) + &
                 fac10(lay) * absa(ind0+1,ig) + &
                 fac01(lay) * absa(ind1,ig) + &
                 fac11(lay) * absa(ind1+1,ig)) &
                 + tauself + taufor &
                 + tauo2
            fracs(lay,ngs10+ig) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(11) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(11) + 1
         indf = indfor(lay)
         indm = indminor(lay)
         scaleo2 = colo2(lay)*scaleminor(lay)
         do ig = 1, ng11
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            tauo2 =  scaleo2 * (kb_mo2(indm,ig) + minorfrac(lay) * &
                 (kb_mo2(indm+1,ig) - kb_mo2(indm,ig)))
            taug(lay,ngs10+ig) = colh2o(lay) * &
                 (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) + &
                 fac11(lay) * absb(ind1+1,ig))  &
                 + taufor &
                 + tauo2
            fracs(lay,ngs10+ig) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb11

!----------------------------------------------------------------------------
      subroutine taugb12
!----------------------------------------------------------------------------
!
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng12, ngs11
      use rrlw_ref, only : chi_mls
      use rrlw_kg12, only : fracrefa, absa, &
                            selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor
      real(kind=r8) :: refrat_planck_a
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb12_codon(nlayers_c, laytrop_c, ng12_c, ngs11_c, &
              nspa12_c, oneminus_c, colh2o_p, colco2_p, rat_h2oco2_p, &
              rat_h2oco2_1_p, jp_p, jt_p, jt1_p, indself_p, indfor_p, &
              fac00_p, fac01_p, fac10_p, fac11_p, selffac_p, selffrac_p, &
              forfac_p, forfrac_p, chi_mls_p, fracrefa_p, absa_p, selfref_p, &
              forref_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb12_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng12_c, ngs11_c
            integer(c_int64_t), value :: nspa12_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: colh2o_p, colco2_p, rat_h2oco2_p, rat_h2oco2_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: chi_mls_p, fracrefa_p, absa_p, selfref_p, forref_p
            type(c_ptr), value :: fracs_p, taug_p
         end subroutine rrtmg_lw_taugb12_codon
      end interface


! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

! P =   174.164 mb 
      refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum adn foreign continuum is interpolated 
! (in temperature) separately.  

      call taugb1_2_6_8_12_16_lw_select_impl()
      if (.not. use_native_taugb1_2_6_8_12_16_lw_impl) then
         call taugb1_2_6_8_12_16_lw_log_entered(12)
         call rrtmg_lw_taugb12_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng12, c_int64_t), &
              int(ngs11, c_int64_t), int(nspa(12), c_int64_t), real(oneminus, c_double), &
              transfer(loc(colh2o(1)), c_null_ptr), transfer(loc(colco2(1)), c_null_ptr), &
              transfer(loc(rat_h2oco2(1)), c_null_ptr), transfer(loc(rat_h2oco2_1(1)), c_null_ptr), &
              c_loc(jp64(1)), c_loc(jt64(1)), c_loc(jt164(1)), c_loc(indself64(1)), &
              c_loc(indfor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(chi_mls(1,1)), c_null_ptr), &
              transfer(loc(fracrefa(1,1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(selfref(1,1)), c_null_ptr), transfer(loc(forref(1,1)), c_null_ptr), &
              transfer(loc(fracs(1,1)), c_null_ptr), transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(12) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(12) + js1
         inds = indself(lay)
         indf = indfor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng12
            tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs11+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracs(lay,ngs11+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         do ig = 1, ng12
            taug(lay,ngs11+ig) = 0.0_r8
            fracs(lay,ngs11+ig) = 0.0_r8
         enddo
      enddo

      end subroutine taugb12

!----------------------------------------------------------------------------
      subroutine taugb13
!----------------------------------------------------------------------------
!
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng13, ngs12
      use rrlw_ref, only : chi_mls
      use rrlw_kg13, only : fracrefa, fracrefb, absa, &
                            ka_mco2, ka_mco, kb_mo3, selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmco2, jmco, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real(kind=r8) :: speccomb_mco, specparm_mco, specmult_mco, fmco
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, co2m1, co2m2, absco2 
      real(kind=r8) :: com1, com2, absco, abso3
      real(kind=r8) :: chi_co2, ratco2, adjfac, adjcolco2
      real(kind=r8) :: refrat_planck_a, refrat_m_a, refrat_m_a3
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb13_codon(nlayers_c, laytrop_c, ng13_c, &
              ngs12_c, nspa13_c, oneminus_c, colh2o_p, coln2o_p, colco2_p, &
              colco_p, colo3_p, coldry_p, rat_h2on2o_p, rat_h2on2o_1_p, &
              jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p, fac00_p, &
              fac01_p, fac10_p, fac11_p, selffac_p, selffrac_p, forfac_p, &
              forfrac_p, minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p, &
              absa_p, ka_mco2_p, ka_mco_p, kb_mo3_p, selfref_p, forref_p, &
              fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb13_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng13_c, ngs12_c
            integer(c_int64_t), value :: nspa13_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: colh2o_p, coln2o_p, colco2_p, colco_p, colo3_p, coldry_p
            type(c_ptr), value :: rat_h2on2o_p, rat_h2on2o_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: minorfrac_p, chi_mls_p, fracrefa_p, fracrefb_p
            type(c_ptr), value :: absa_p, ka_mco2_p, ka_mco_p, kb_mo3_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb13_codon
      end interface


! Minor gas mapping levels :
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - co, p = 706 mb, t = 278.94 k
!     upper - o3, p = 95.5835 mb, t = 215.7 k

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower/upper atmosphere.

! P = 473.420 mb (Level 5)
      refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)

! P = 1053. (Level 1)
      refrat_m_a = chi_mls(1,1)/chi_mls(4,1)

! P = 706. (Level 3)
      refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

      call taugb3_7_9_13_15_lw_select_impl()
      if (.not. use_native_taugb3_7_9_13_15_lw_impl) then
         call taugb3_7_9_13_15_lw_log_entered(13)
         call rrtmg_lw_taugb13_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng13, c_int64_t), &
              int(ngs12, c_int64_t), int(nspa(13), c_int64_t), real(oneminus, c_double), &
              transfer(loc(colh2o(1)), c_null_ptr), transfer(loc(coln2o(1)), c_null_ptr), &
              transfer(loc(colco2(1)), c_null_ptr), transfer(loc(colco(1)), c_null_ptr), &
              transfer(loc(colo3(1)), c_null_ptr), transfer(loc(coldry(1)), c_null_ptr), &
              transfer(loc(rat_h2on2o(1)), c_null_ptr), transfer(loc(rat_h2on2o_1(1)), c_null_ptr), &
              c_loc(jp64(1)), c_loc(jt64(1)), c_loc(jt164(1)), c_loc(indself64(1)), &
              c_loc(indfor64(1)), c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(minorfrac(1)), c_null_ptr), &
              transfer(loc(chi_mls(1,1)), c_null_ptr), transfer(loc(fracrefa(1,1)), c_null_ptr), &
              transfer(loc(fracrefb(1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(ka_mco2(1,1,1)), c_null_ptr), transfer(loc(ka_mco(1,1,1)), c_null_ptr), &
              transfer(loc(kb_mo3(1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2on2o(lay)*coln2o(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(lay) + rat_h2on2o_1(lay)*coln2o(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mco2 = colh2o(lay) + refrat_m_a*coln2o(lay)
         specparm_mco2 = colh2o(lay)/speccomb_mco2
         if (specparm_mco2 .ge. oneminus) specparm_mco2 = oneminus
         specmult_mco2 = 8._r8*specparm_mco2
         jmco2 = 1 + int(specmult_mco2)
         fmco2 = mod(specmult_mco2,1.0_r8)

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor 
!  to obtain the proper contribution.
         chi_co2 = colco2(lay)/(coldry(lay))
         ratco2 = 1.e20_r8*chi_co2/3.55e-4_r8
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.68_r8
            adjcolco2 = adjfac*3.55e-4*coldry(lay)*1.e-20_r8
         else
            adjcolco2 = colco2(lay)
         endif

         speccomb_mco = colh2o(lay) + refrat_m_a3*coln2o(lay)
         specparm_mco = colh2o(lay)/speccomb_mco
         if (specparm_mco .ge. oneminus) specparm_mco = oneminus
         specmult_mco = 8._r8*specparm_mco
         jmco = 1 + int(specmult_mco)
         fmco = mod(specmult_mco,1.0_r8)

         speccomb_planck = colh2o(lay)+refrat_planck_a*coln2o(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(13) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(13) + js1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng13
            tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1)
            com1 = ka_mco(jmco,indm,ig) + fmco * &
                 (ka_mco(jmco+1,indm,ig) - ka_mco(jmco,indm,ig))
            com2 = ka_mco(jmco,indm+1,ig) + fmco * &
                 (ka_mco(jmco+1,indm+1,ig) - ka_mco(jmco,indm+1,ig))
            absco = com1 + minorfrac(lay) * (com2 - com1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs12+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colco(lay)*absco
            fracs(lay,ngs12+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         indm = indminor(lay)
         do ig = 1, ng13
            abso3 = kb_mo3(indm,ig) + minorfrac(lay) * &
                 (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
            taug(lay,ngs12+ig) = colo3(lay)*abso3
            fracs(lay,ngs12+ig) =  fracrefb(ig)
         enddo
      enddo

      end subroutine taugb13

!----------------------------------------------------------------------------
      subroutine taugb14
!----------------------------------------------------------------------------
!
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng14, ngs13
      use rrlw_kg14, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, ig
      real(kind=r8) :: tauself, taufor

      interface
         subroutine rrtmg_lw_taugb14_codon(nlayers_c, laytrop_c, ng14_c, ngs13_c, &
              nspa14_c, nspb14_c, colco2_p, jp_p, jt_p, jt1_p, indself_p, &
              indfor_p, fac00_p, fac01_p, fac10_p, fac11_p, selffac_p, &
              selffrac_p, forfac_p, forfrac_p, fracrefa_p, fracrefb_p, absa_p, &
              absb_p, selfref_p, forref_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb14_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng14_c, ngs13_c
            integer(c_int64_t), value :: nspa14_c, nspb14_c
            type(c_ptr), value :: colco2_p, jp_p, jt_p, jt1_p, indself_p, indfor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p, selffac_p, selffrac_p
            type(c_ptr), value :: forfac_p, forfrac_p, fracrefa_p, fracrefb_p, absa_p, absb_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb14_codon
      end interface

! Compute the optical depth by interpolating in ln(pressure) and 
! temperature.  Below laytrop, the water vapor self-continuum 
! and foreign continuum is interpolated (in temperature) separately.  

      call taugb10_11_14_lw_select_impl()
      if (.not. use_native_taugb10_11_14_lw_impl) then
         call taugb10_11_14_lw_log_entered(14)
         call rrtmg_lw_taugb14_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng14, c_int64_t), &
              int(ngs13, c_int64_t), int(nspa(14), c_int64_t), int(nspb(14), c_int64_t), &
              transfer(loc(colco2(1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              transfer(loc(fac00(1)), c_null_ptr), transfer(loc(fac01(1)), c_null_ptr), &
              transfer(loc(fac10(1)), c_null_ptr), transfer(loc(fac11(1)), c_null_ptr), &
              transfer(loc(selffac(1)), c_null_ptr), transfer(loc(selffrac(1)), c_null_ptr), &
              transfer(loc(forfac(1)), c_null_ptr), transfer(loc(forfrac(1)), c_null_ptr), &
              transfer(loc(fracrefa(1)), c_null_ptr), transfer(loc(fracrefb(1)), c_null_ptr), &
              transfer(loc(absa(1,1)), c_null_ptr), transfer(loc(absb(1,1)), c_null_ptr), &
              transfer(loc(selfref(1,1)), c_null_ptr), transfer(loc(forref(1,1)), c_null_ptr), &
              transfer(loc(fracs(1,1)), c_null_ptr), transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop
         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(14) + 1
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(14) + 1
         inds = indself(lay)
         indf = indfor(lay)
         do ig = 1, ng14
            tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            taug(lay,ngs13+ig) = colco2(lay) * &
                 (fac00(lay) * absa(ind0,ig) + &
                 fac10(lay) * absa(ind0+1,ig) + &
                 fac01(lay) * absa(ind1,ig) + &
                 fac11(lay) * absa(ind1+1,ig)) &
                 + tauself + taufor
            fracs(lay,ngs13+ig) = fracrefa(ig)
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(14) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(14) + 1
         do ig = 1, ng14
            taug(lay,ngs13+ig) = colco2(lay) * &
                 (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) + &
                 fac11(lay) * absb(ind1+1,ig))
            fracs(lay,ngs13+ig) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb14

!----------------------------------------------------------------------------
      subroutine taugb15
!----------------------------------------------------------------------------
!
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng15, ngs14
      use rrlw_ref, only : chi_mls
      use rrlw_kg15, only : fracrefa, absa, &
                            ka_mn2, selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmn2, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mn2, specparm_mn2, specmult_mn2, fmn2
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: scalen2, tauself, taufor, n2m1, n2m2, taun2 
      real(kind=r8) :: refrat_planck_a, refrat_m_a
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb15_codon(nlayers_c, laytrop_c, ng15_c, &
              ngs14_c, nspa15_c, oneminus_c, coln2o_p, colco2_p, colbrd_p, &
              rat_n2oco2_p, rat_n2oco2_1_p, jp_p, jt_p, jt1_p, indself_p, &
              indfor_p, indminor_p, fac00_p, fac01_p, fac10_p, fac11_p, &
              selffac_p, selffrac_p, forfac_p, forfrac_p, scaleminor_p, &
              minorfrac_p, chi_mls_p, fracrefa_p, absa_p, ka_mn2_p, &
              selfref_p, forref_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb15_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng15_c, ngs14_c
            integer(c_int64_t), value :: nspa15_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: coln2o_p, colco2_p, colbrd_p
            type(c_ptr), value :: rat_n2oco2_p, rat_n2oco2_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p, indminor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: scaleminor_p, minorfrac_p, chi_mls_p, fracrefa_p
            type(c_ptr), value :: absa_p, ka_mn2_p, selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb15_codon
      end interface


! Minor gas mapping level : 
!     Lower - Nitrogen Continuum, P = 1053., T = 294.

! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.
! P = 1053. mb (Level 1)
      refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)

! P = 1053.
      refrat_m_a = chi_mls(4,1)/chi_mls(2,1)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

      call taugb3_7_9_13_15_lw_select_impl()
      if (.not. use_native_taugb3_7_9_13_15_lw_impl) then
         call taugb3_7_9_13_15_lw_log_entered(15)
         call rrtmg_lw_taugb15_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng15, c_int64_t), &
              int(ngs14, c_int64_t), int(nspa(15), c_int64_t), real(oneminus, c_double), &
              transfer(loc(coln2o(1)), c_null_ptr), transfer(loc(colco2(1)), c_null_ptr), &
              transfer(loc(colbrd(1)), c_null_ptr), transfer(loc(rat_n2oco2(1)), c_null_ptr), &
              transfer(loc(rat_n2oco2_1(1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              c_loc(indminor64(1)), transfer(loc(fac00(1)), c_null_ptr), &
              transfer(loc(fac01(1)), c_null_ptr), transfer(loc(fac10(1)), c_null_ptr), &
              transfer(loc(fac11(1)), c_null_ptr), transfer(loc(selffac(1)), c_null_ptr), &
              transfer(loc(selffrac(1)), c_null_ptr), transfer(loc(forfac(1)), c_null_ptr), &
              transfer(loc(forfrac(1)), c_null_ptr), transfer(loc(scaleminor(1)), c_null_ptr), &
              transfer(loc(minorfrac(1)), c_null_ptr), transfer(loc(chi_mls(1,1)), c_null_ptr), &
              transfer(loc(fracrefa(1,1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(ka_mn2(1,1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = coln2o(lay) + rat_n2oco2(lay)*colco2(lay)
         specparm = coln2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = coln2o(lay) + rat_n2oco2_1(lay)*colco2(lay)
         specparm1 = coln2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mn2 = coln2o(lay) + refrat_m_a*colco2(lay)
         specparm_mn2 = coln2o(lay)/speccomb_mn2
         if (specparm_mn2 .ge. oneminus) specparm_mn2 = oneminus
         specmult_mn2 = 8._r8*specparm_mn2
         jmn2 = 1 + int(specmult_mn2)
         fmn2 = mod(specmult_mn2,1.0_r8)

         speccomb_planck = coln2o(lay)+refrat_planck_a*colco2(lay)
         specparm_planck = coln2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(15) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(15) + js1
         inds = indself(lay)
         indf = indfor(lay)
         indm = indminor(lay)
         
         scalen2 = colbrd(lay)*scaleminor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif
         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng15
            tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 
            n2m1 = ka_mn2(jmn2,indm,ig) + fmn2 * &
                 (ka_mn2(jmn2+1,indm,ig) - ka_mn2(jmn2,indm,ig))
            n2m2 = ka_mn2(jmn2,indm+1,ig) + fmn2 * &
                 (ka_mn2(jmn2+1,indm+1,ig) - ka_mn2(jmn2,indm+1,ig))
            taun2 = scalen2 * (n2m1 + minorfrac(lay) * (n2m2 - n2m1))

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif 

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs14+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + taun2
            fracs(lay,ngs14+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         do ig = 1, ng15
            taug(lay,ngs14+ig) = 0.0_r8
            fracs(lay,ngs14+ig) = 0.0_r8
         enddo
      enddo

      end subroutine taugb15

!----------------------------------------------------------------------------
      subroutine taugb16
!----------------------------------------------------------------------------
!
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!----------------------------------------------------------------------------

! ------- Modules -------

      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use parrrtm, only : ng16, ngs15
      use rrlw_ref, only : chi_mls
      use rrlw_kg16, only : fracrefa, fracrefb, absa, absb, &
                            selfref, forref

! ------- Declarations -------

! Local 
      integer :: lay, ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor
      real(kind=r8) :: refrat_planck_a
      real(kind=r8) :: tau_major, tau_major1

      interface
         subroutine rrtmg_lw_taugb16_codon(nlayers_c, laytrop_c, ng16_c, ngs15_c, &
              nspa16_c, nspb16_c, oneminus_c, colh2o_p, colch4_p, rat_h2och4_p, &
              rat_h2och4_1_p, jp_p, jt_p, jt1_p, indself_p, indfor_p, fac00_p, &
              fac01_p, fac10_p, fac11_p, selffac_p, selffrac_p, forfac_p, &
              forfrac_p, chi_mls_p, fracrefa_p, fracrefb_p, absa_p, absb_p, &
              selfref_p, forref_p, fracs_p, taug_p) bind(c, name="rrtmg_lw_taugb16_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlayers_c, laytrop_c, ng16_c, ngs15_c
            integer(c_int64_t), value :: nspa16_c, nspb16_c
            real(c_double), value :: oneminus_c
            type(c_ptr), value :: colh2o_p, colch4_p, rat_h2och4_p, rat_h2och4_1_p
            type(c_ptr), value :: jp_p, jt_p, jt1_p, indself_p, indfor_p
            type(c_ptr), value :: fac00_p, fac01_p, fac10_p, fac11_p
            type(c_ptr), value :: selffac_p, selffrac_p, forfac_p, forfrac_p
            type(c_ptr), value :: chi_mls_p, fracrefa_p, fracrefb_p, absa_p, absb_p
            type(c_ptr), value :: selfref_p, forref_p, fracs_p, taug_p
         end subroutine rrtmg_lw_taugb16_codon
      end interface


! Calculate reference ratio to be used in calculation of Planck
! fraction in lower atmosphere.

! P = 387. mb (Level 6)
      refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)

! Compute the optical depth by interpolating in ln(pressure), 
! temperature,and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated 
! (in temperature) separately.  

      call taugb1_2_6_8_12_16_lw_select_impl()
      if (.not. use_native_taugb1_2_6_8_12_16_lw_impl) then
         call taugb1_2_6_8_12_16_lw_log_entered(16)
         call rrtmg_lw_taugb16_codon( &
              int(nlayers, c_int64_t), int(laytrop, c_int64_t), int(ng16, c_int64_t), &
              int(ngs15, c_int64_t), int(nspa(16), c_int64_t), int(nspb(16), c_int64_t), &
              real(oneminus, c_double), transfer(loc(colh2o(1)), c_null_ptr), &
              transfer(loc(colch4(1)), c_null_ptr), transfer(loc(rat_h2och4(1)), c_null_ptr), &
              transfer(loc(rat_h2och4_1(1)), c_null_ptr), c_loc(jp64(1)), c_loc(jt64(1)), &
              c_loc(jt164(1)), c_loc(indself64(1)), c_loc(indfor64(1)), &
              transfer(loc(fac00(1)), c_null_ptr), transfer(loc(fac01(1)), c_null_ptr), &
              transfer(loc(fac10(1)), c_null_ptr), transfer(loc(fac11(1)), c_null_ptr), &
              transfer(loc(selffac(1)), c_null_ptr), transfer(loc(selffrac(1)), c_null_ptr), &
              transfer(loc(forfac(1)), c_null_ptr), transfer(loc(forfrac(1)), c_null_ptr), &
              transfer(loc(chi_mls(1,1)), c_null_ptr), transfer(loc(fracrefa(1,1)), c_null_ptr), &
              transfer(loc(fracrefb(1)), c_null_ptr), transfer(loc(absa(1,1)), c_null_ptr), &
              transfer(loc(absb(1,1)), c_null_ptr), transfer(loc(selfref(1,1)), c_null_ptr), &
              transfer(loc(forref(1,1)), c_null_ptr), transfer(loc(fracs(1,1)), c_null_ptr), &
              transfer(loc(taug(1,1)), c_null_ptr) &
         )
         return
      endif

! Lower atmosphere loop
      do lay = 1, laytrop

         speccomb = colh2o(lay) + rat_h2och4(lay)*colch4(lay)
         specparm = colh2o(lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(lay) + rat_h2och4_1(lay)*colch4(lay)
         specparm1 = colh2o(lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay)
         specparm_planck = colh2o(lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(16) + js
         ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(16) + js1
         inds = indself(lay)
         indf = indfor(lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(lay)
            fac100 = fk1*fac00(lay)
            fac200 = fk2*fac00(lay)
            fac010 = fk0*fac10(lay)
            fac110 = fk1*fac10(lay)
            fac210 = fk2*fac10(lay)
         else
            fac000 = (1._r8 - fs) * fac00(lay)
            fac010 = (1._r8 - fs) * fac10(lay)
            fac100 = fs * fac00(lay)
            fac110 = fs * fac10(lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1 
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(lay)
            fac101 = fk1*fac01(lay)
            fac201 = fk2*fac01(lay)
            fac011 = fk0*fac11(lay)
            fac111 = fk1*fac11(lay)
            fac211 = fk2*fac11(lay)
         else
            fac001 = (1._r8 - fs1) * fac01(lay)
            fac011 = (1._r8 - fs1) * fac11(lay)
            fac101 = fs1 * fac01(lay)
            fac111 = fs1 * fac11(lay)
         endif

         do ig = 1, ng16
            tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                 (forref(indf+1,ig) - forref(indf,ig))) 

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac200 * absa(ind0+2,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig) + &
                    fac210 * absa(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa(ind0-1,ig) + &
                    fac100 * absa(ind0,ig) + &
                    fac000 * absa(ind0+1,ig) + &
                    fac210 * absa(ind0+8,ig) + &
                    fac110 * absa(ind0+9,ig) + &
                    fac010 * absa(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa(ind0,ig) + &
                    fac100 * absa(ind0+1,ig) + &
                    fac010 * absa(ind0+9,ig) + &
                    fac110 * absa(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac201 * absa(ind1+2,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig) + &
                    fac211 * absa(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa(ind1-1,ig) + &
                    fac101 * absa(ind1,ig) + &
                    fac001 * absa(ind1+1,ig) + &
                    fac211 * absa(ind1+8,ig) + &
                    fac111 * absa(ind1+9,ig) + &
                    fac011 * absa(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa(ind1,ig) + &
                    fac101 * absa(ind1+1,ig) + &
                    fac011 * absa(ind1+9,ig) + &
                    fac111 * absa(ind1+10,ig))
            endif

            taug(lay,ngs15+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracs(lay,ngs15+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
         enddo
      enddo

! Upper atmosphere loop
      do lay = laytrop+1, nlayers
         ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(16) + 1
         ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(16) + 1
         do ig = 1, ng16
            taug(lay,ngs15+ig) = colch4(lay) * &
                 (fac00(lay) * absb(ind0,ig) + &
                 fac10(lay) * absb(ind0+1,ig) + &
                 fac01(lay) * absb(ind1,ig) + &
                 fac11(lay) * absb(ind1+1,ig))
            fracs(lay,ngs15+ig) = fracrefb(ig)
         enddo
      enddo

      end subroutine taugb16

      end subroutine taumol

! --------------------------------------------------------------------------
      subroutine taugb4_5_lw_select_impl()

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (taugb4_5_lw_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('RRTMG_LW_TAUGB4_5_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_taugb4_5_lw_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_taugb4_5_lw_impl = .false.
      end if

      taugb4_5_lw_impl_selected = .true.

      if (masterproc) then
         if (use_native_taugb4_5_lw_impl) then
            write(iulog,*) 'rrtmg_lw_taugb4_5 implementation = native'
         else
            write(iulog,*) 'rrtmg_lw_taugb4_5 implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine taugb4_5_lw_select_impl

! --------------------------------------------------------------------------
      subroutine taugb4_5_lw_log_entered(band)

      integer, intent(in) :: band
      integer :: idx
      character(len=96) :: proof_line

      select case (band)
      case (4)
         idx = 1
         proof_line = 'rrtmg_lw_taugb4 entered (longwave band 4 optical depth = codon)'
      case (5)
         idx = 2
         proof_line = 'rrtmg_lw_taugb5 entered (longwave band 5 optical depth = codon)'
      case default
         return
      end select

      if (taugb4_5_lw_entered_logged(idx)) return
      taugb4_5_lw_entered_logged(idx) = .true.

      if (masterproc) then
         write(iulog,*) trim(proof_line)
         call flush(iulog)
      end if

      end subroutine taugb4_5_lw_log_entered

! --------------------------------------------------------------------------
      subroutine taugb3_7_9_13_15_lw_select_impl()

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (taugb3_7_9_13_15_lw_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('RRTMG_LW_TAUGB3_7_9_13_15_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_taugb3_7_9_13_15_lw_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_taugb3_7_9_13_15_lw_impl = .false.
      end if

      taugb3_7_9_13_15_lw_impl_selected = .true.

      if (masterproc) then
         if (use_native_taugb3_7_9_13_15_lw_impl) then
            write(iulog,*) 'rrtmg_lw_taugb3_7_9_13_15 implementation = native'
         else
            write(iulog,*) 'rrtmg_lw_taugb3_7_9_13_15 implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine taugb3_7_9_13_15_lw_select_impl

! --------------------------------------------------------------------------
      subroutine taugb3_7_9_13_15_lw_log_entered(band)

      integer, intent(in) :: band
      integer :: idx
      character(len=96) :: proof_line

      select case (band)
      case (3)
         idx = 1
         proof_line = 'rrtmg_lw_taugb3 entered (longwave band 3 optical depth = codon)'
      case (7)
         idx = 2
         proof_line = 'rrtmg_lw_taugb7 entered (longwave band 7 optical depth = codon)'
      case (9)
         idx = 3
         proof_line = 'rrtmg_lw_taugb9 entered (longwave band 9 optical depth = codon)'
      case (13)
         idx = 4
         proof_line = 'rrtmg_lw_taugb13 entered (longwave band 13 optical depth = codon)'
      case (15)
         idx = 5
         proof_line = 'rrtmg_lw_taugb15 entered (longwave band 15 optical depth = codon)'
      case default
         return
      end select

      if (taugb3_7_9_13_15_lw_entered_logged(idx)) return
      taugb3_7_9_13_15_lw_entered_logged(idx) = .true.

      if (masterproc) then
         write(iulog,*) trim(proof_line)
         call flush(iulog)
      end if

      end subroutine taugb3_7_9_13_15_lw_log_entered

! --------------------------------------------------------------------------
      subroutine taugb1_2_6_8_12_16_lw_select_impl()

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (taugb1_2_6_8_12_16_lw_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('RRTMG_LW_TAUGB1_2_6_8_12_16_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_taugb1_2_6_8_12_16_lw_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_taugb1_2_6_8_12_16_lw_impl = .false.
      end if

      taugb1_2_6_8_12_16_lw_impl_selected = .true.

      if (masterproc) then
         if (use_native_taugb1_2_6_8_12_16_lw_impl) then
            write(iulog,*) 'rrtmg_lw_taugb1_2_6_8_12_16 implementation = native'
         else
            write(iulog,*) 'rrtmg_lw_taugb1_2_6_8_12_16 implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine taugb1_2_6_8_12_16_lw_select_impl

! --------------------------------------------------------------------------
      subroutine taugb1_2_6_8_12_16_lw_log_entered(band)

      integer, intent(in) :: band
      integer :: idx
      character(len=96) :: proof_line

      select case (band)
      case (1)
         idx = 1
         proof_line = 'rrtmg_lw_taugb1 entered (longwave band 1 optical depth = codon)'
      case (2)
         idx = 2
         proof_line = 'rrtmg_lw_taugb2 entered (longwave band 2 optical depth = codon)'
      case (6)
         idx = 3
         proof_line = 'rrtmg_lw_taugb6 entered (longwave band 6 optical depth = codon)'
      case (8)
         idx = 4
         proof_line = 'rrtmg_lw_taugb8 entered (longwave band 8 optical depth = codon)'
      case (12)
         idx = 5
         proof_line = 'rrtmg_lw_taugb12 entered (longwave band 12 optical depth = codon)'
      case (16)
         idx = 6
         proof_line = 'rrtmg_lw_taugb16 entered (longwave band 16 optical depth = codon)'
      case default
         return
      end select

      if (taugb1_2_6_8_12_16_lw_entered_logged(idx)) return
      taugb1_2_6_8_12_16_lw_entered_logged(idx) = .true.

      if (masterproc) then
         write(iulog,*) trim(proof_line)
         call flush(iulog)
      end if

      end subroutine taugb1_2_6_8_12_16_lw_log_entered

! --------------------------------------------------------------------------
      subroutine taugb10_11_14_lw_select_impl()

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (taugb10_11_14_lw_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('RRTMG_LW_TAUGB10_11_14_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_taugb10_11_14_lw_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_taugb10_11_14_lw_impl = .false.
      end if

      taugb10_11_14_lw_impl_selected = .true.

      if (masterproc) then
         if (use_native_taugb10_11_14_lw_impl) then
            write(iulog,*) 'rrtmg_lw_taugb10_11_14 implementation = native'
         else
            write(iulog,*) 'rrtmg_lw_taugb10_11_14 implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine taugb10_11_14_lw_select_impl

! --------------------------------------------------------------------------
      subroutine taugb10_11_14_lw_log_entered(band)

      integer, intent(in) :: band
      integer :: idx
      character(len=96) :: proof_line

      select case (band)
      case (10)
         idx = 1
         proof_line = 'rrtmg_lw_taugb10 entered (longwave band 10 optical depth = codon)'
      case (11)
         idx = 2
         proof_line = 'rrtmg_lw_taugb11 entered (longwave band 11 optical depth = codon)'
      case (14)
         idx = 3
         proof_line = 'rrtmg_lw_taugb14 entered (longwave band 14 optical depth = codon)'
      case default
         return
      end select

      if (taugb10_11_14_lw_entered_logged(idx)) return
      taugb10_11_14_lw_entered_logged(idx) = .true.

      if (masterproc) then
         write(iulog,*) trim(proof_line)
         call flush(iulog)
      end if

      end subroutine taugb10_11_14_lw_log_entered

      subroutine taumol_lw_log_impl()

      call taugb1_2_6_8_12_16_lw_select_impl()
      call taugb3_7_9_13_15_lw_select_impl()
      call taugb4_5_lw_select_impl()
      call taugb10_11_14_lw_select_impl()
      if (use_native_taugb1_2_6_8_12_16_lw_impl .or. use_native_taugb3_7_9_13_15_lw_impl .or. &
          use_native_taugb4_5_lw_impl .or. use_native_taugb10_11_14_lw_impl) return

      if (taumol_lw_impl_logged) return
      taumol_lw_impl_logged = .true.

      if (masterproc) then
         write(iulog,*) 'taumol implementation = codon'
         call flush(iulog)
      end if

      end subroutine taumol_lw_log_impl

      end module rrtmg_lw_taumol
