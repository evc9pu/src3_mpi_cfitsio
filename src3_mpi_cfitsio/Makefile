#FC = mpifrtpx 
#FC = mpif90 -DMPI -ffree-line-length-none
FC = mpif90 -DMPI -ffree-line-length-none -fallow-argument-mismatch
#FC = mpif90 -DMPI -ffree-line-length-none -fno-strict-aliasing
#FC = mpifort -DMPI -ffree-line-length-none -fallow-argument-mismatch
#FC = mpif90
#FC = mpif90 -DMPI
#FC = gfortran
#CFITSIO = -mcmodel=medium
CFITSIO = -L/usr/local/lib/ -lcfitsio
LDFLAGS = $(CFITSIO)

.SUFFIXES:
.SUFFIXES: .o .f .f90

OBJSE	= 	 \
		math_binning.o \
		type_yale_sparse_matrix.o type_smart_image.o type_nested_image.o \
		lib_array.o constants.o \
		extract_nested_image.o

OBJSL	= 	 \
		math_binning.o \
		type_yale_sparse_matrix.o type_smart_image.o type_nested_image.o \
		nested_image_slice.o

OBJSB	=	constants.o \
		messages.o \
		lib_array.o \
		atmos_interp_mod.o \
		lib_cfitsio.o\
		closest_wall.o \
		math_binning.o type_yale_sparse_matrix.o \
		ttsre_mpi_mod.o \
		type_smart_image.o type_nested_image_nofits.o \
		grid_mod_5.o tts_mod_5.o dust_mod_5.o stokes_mod.o filt_mod.o \
		taunum_mod.o out_mod.o opacin_mod_2.o log_mod.o tab_mod.o \
		ecuyer_cote_mod.o \
		random_standard_uniform_mod.o \
		random.o tabl_mod.o spot_mod.o isrf_mod.o draine_mod.o \
		user_set_generator.o \
		tempmap_mod.o \
		lineint_mod.o \
		depint_mod.o \
		peeloff_3d_2_test.o \
		tauint_3d_2_fuv.o \
		tauint_peel.o \
		output_FITS_2_test.o \
		alphadisk_new.o \
		check.o \
		configuration.o \
		cop.o \
		densdiskg_new1_model12_13_evo.o \
		densenv.o \
		densenvpw.o \
		densenv_read_5.o \
		densdiskwind_new1.o \
		depletion_co.o \
		depmap.o \
		deppeel.o \
		diskprofile_model11.o \
		draine_cdf.o \
		emit_sg.o \
		emit_gasacc.o \
		flux_test.o \
		dust_5.o \
		dustmat.o \
		dusttemp_5.o \
		envset.o \
		errmsg.o \
		fractal_sph.o \
		filt_nu.o \
		findangle3.o \
		findopacid_5.o \
		fuv_direct.o \
		gridset_5_model13new.o \
		initp_common.o \
		initp.o \
		initp_ps.o \
		initpspot.o \
		initpacc_4.o \
		initpabs.o \
		initpout.o \
		int_mean.o \
		is_fuv.o \
		isnan.o \
		isrf_cdf.o \
		isotrp.o \
		linemap_old.o \
		linepeel.o \
		locate3.o \
		lucydustfreq_2_test.o \
		mathsubs.o \
		newdisktherm_5_fuv.o \
		nr_2.o \
		phiface.o \
		plancknu.o \
		photosphere.o \
		photospherepeel.o \
		propagate_5_fuv.o \
		radface.o \
		reapar_4_test1.o \
		reemit.o \
		rep_isrf.o \
		rep_draine_new.o \
		samptabl.o \
		select_dust_2.o \
		setdiffus_5.o \
		setup.o \
		setup_wave.o \
        spotset.o \
		splines.o \
		starbarb_2.o \
		stokes.o \
		stokespeel.o \
		subs_g77.o \
        Rdist.o \
		testgrid.o \
		tempmap.o \
		tempeel.o \
		tfinal_5_test1.o \
		thetaface.o \
		wrimsg.o \
		zerod.o 

OBJSB1  =	emit_common_2.o \
		initarr_2_new.o \
		newtts.o 

OBJSB2  =	emit_common_final.o \
		initarr_final.o \
		newtts_final.o \


ttsre:		$(OBJSB) $(OBJSB1)
		$(FC) $(OBJSB) $(OBJSB1) $(LDFLAGS) -o p_ttsre_model13new

ttsre_f:	$(OBJSB) $(OBJSB2)
		$(FC) $(OBJSB) $(OBJSB2) $(LDFLAGS) -o p_ttsre_model13new_final


nested_image_slice:	$(OBJSL)
			$(FC) $(OBJSL) $(LDFLAGS) -o nested_image_slice

extract_nested_image:	$(OBJSE)
			$(FC) $(OBJSE) $(LDFLAGS) -o extract_nested_image


.f90.o:
		$(FC) -c $<

.f.o:
		$(FC) -c $<

tarfile:;	tar cvf ttsre.tar *.f90 Make* 

clean:;		/bin/rm -f *.o */*.o */*.mod *genmod.f90 *.mod

