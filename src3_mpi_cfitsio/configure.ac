AC_INIT(release.2009,version-1.0)

# Check for fortran compiler

AC_PROG_FC([ifort f95 gfortran xlf g95])

# Extra options for each compiler

case $FC in
f95)
	AC_SUBST(extra,'-maxcontin=50 -f2003 -colour -ieee=full');;
ifort)
	AC_SUBST(extra,'-heap-arrays');;
g95)
	AC_SUBST(extra,'-ffree-line-length-huge');;
gfortran)
	AC_SUBST(extra,'-ffree-line-length-none');;
*)
	AC_SUBST(extra,'');;
esac

# Check for FITSIO if requested

AC_ARG_ENABLE(fits,[  --enable-fits           Replaces ASCII output with FITS output])

if test "$enable_fits" == "yes"
then

	AC_MSG_NOTICE([Using FITS output])

	echo 'checking for cfitsio... \c'

	AC_ARG_WITH(cfitsio,
		[  --with-cfitsio=<path>   path containing libcfitsio.a],
		[CFITSIO_lib_check="$with_cfitsio"],
		[CFITSIO_lib_check="/usr/local/lib/ /opt/local/lib/ /usr/lib/ /sw/lib/"])

	CFITSIO_libdir=
	for dir in $CFITSIO_lib_check
	do
		if test -d "$dir" && ( test -f "$dir/libCFITSIO.a" )
		then
			CFITSIO_libdir=$dir
			break
		fi
	done

	if test -z "$CFITSIO_libdir"
	then
		echo 'no'
		AC_MSG_ERROR([If you have libcfitsio.a installed, please use the --with-cfitsio=path option to specify its parent directory])
	fi
		
	echo "$CFITSIO_libdir"/libcfitsio.a
	
	AC_SUBST(output,'output_FITS.o')
	AC_SUBST(nested_image,'type_nested_image.o')
	AC_SUBST(CFITSIO_libdir,"-L$CFITSIO_libdir -lcfitsio")
	AC_SUBST(fits,'lib_cfitsio.o')
	
else
	AC_MSG_NOTICE([Using ASCII output (specify --enable-fits to enable FITS file output)])
	AC_SUBST(output,'output.o')
	AC_SUBST(nested_image,'type_nested_image_nofits.o')
	AC_SUBST(CFITSIO_libdir,'')
	AC_SUBST(fits,'')
fi

# Check for debug mode

AC_ARG_ENABLE(debug,[  --enable-debug          Enables the full debugging flags for the particular compiler])

if test "$enable_debug" == "yes"
then
	
	AC_MSG_NOTICE([Using flags for debugging])
	
	case $FC in
	f95)
		AC_SUBST(debug,'-C=all -gline');;
	ifort)
		AC_SUBST(debug,'-check all -warn all -warn nodec -warn interfaces -traceback -g');;
	gfortran)
		AC_SUBST(debug,'-Wall -fbounds-check -g');;
	g95)
		AC_SUBST(debug,'-Wall -fbounds-check -g -ftrace=full');;
	*)
		AC_MSG_NOTICE(["-> Debug flags not implemented for this compiler"])
		AC_SUBST(debug,'');;
	esac
	
else

	AC_SUBST(debug,'')

fi

# Output Makefile

AC_OUTPUT(Makefile)
