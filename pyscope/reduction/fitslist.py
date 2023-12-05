#!/usr/bin/env python

import glob
import logging

import click
import prettytable
from astropy import coordinates as coord
from astropy import time
from astropy.io import fits
import numpy as np

logger = logging.getLogger(__name__)

def round_exptime(exptime):
    """If the exposure time is greater than 10, 
    round to the nearest 1 second. Otherwise,
    if greater than 1, round to the nearest 0.1 second, 
    otherwise round to 0.001 seconds.
    """
    if exptime > 10:
        return round(exptime)
    elif exptime > 1:
        return round(exptime, 1)
    else:
        return round(exptime, 3)


@click.command(
    epilog="""Check out the documentation at
                https://pyscope.readthedocs.io/ for more
                information."""
)
@click.option("-d", "--date", "check_date", default="", help="Date [default all].")
@click.option("-f", "--filt", "check_filt", default="", help="Filter name [default all].")
@click.option("-r", "--readout", "check_readout", default="", help="Readout mode [default all].")
@click.option("-b", "--binning", "check_binning", default="", help="Binning [default all].")
@click.option(
    "-e",
    "--exptime",
    "check_exptime",
    default="",
    help=f"""Approximate exposure time [default all].
                Note that an error of up to 1% is permitted to allow for imprecisions
                in the camera.""",
)
@click.option("-t", "--target", "check_target", default="", help="Target name [default all].")
@click.option(
    "-v", "--verbose", count=True, type=click.IntRange(0, 1), help="Verbose output."
)
@click.argument("fnames", nargs=-1, type=click.Path(exists=True, dir_okay=False))
@click.version_option()

def fitslist_cli(
    check_date,
    check_filt,
    check_readout,
    check_binning,
    check_exptime,
    check_target,
    verbose,
    fnames,
):
    """List FITS files and their properties."""

    # Set up logging
    logger.setLevel(int(10 * (1 - verbose)))
    logger.debug(
        f"""filt={check_filt}, readout={check_readout}, binning={check_binning}, exptime={check_exptime},
                target={check_target}, verbose={verbose}, fnames={fnames}"""
    )

    # Get list of files
    ftsfiles = []
    if len(fnames) == 0:
        logger.debug("No arguments passed. Using glob to find files.")
        for ext in (".fts", ".fits", ".fit"):
            ftsfiles.extend(glob.glob(f"./*{ext}"))
    else:
        ftsfiles = fnames
    logger.debug(f"fnames={fnames}")
    logger.debug(f"Found {len(fnames)} files.")

    print_rows = []
    for ftsfile in ftsfiles:
        try:
            with fits.open(ftsfile) as hdu:
                header = hdu[0].header
                # data = hdu[0].data # Currently unused
        except:
            logger.warning(f"Could not open {ftsfile}.")
            continue

        # Get properties
        date = time.Time(header["DATE-OBS"], format="fits", scale="utc")
        if date.strftime("%Y-%m-%d") not in check_date.split(",") or check_date == "":
            logger.debug(
                f"Date {date.strftime('%Y-%m-%d')} not in {check_date}. Skipping {ftsfile}."
            )
        
        # Filter
        try:
            filt = header["FILTER"]
        except KeyError:
            try:
                filt = header["FILT"]
            except KeyError:
                filt = ""
                click.echo("Can't find filter.")
        if filt not in check_filt.split(",") or check_filt == "":
            logger.debug(f"Filter {filt} not in {check_filt}. Skipping {ftsfile}.")

        # Readout mode
        try:
            readout_mode = header["READOUTM"]
        except KeyError:
            try:
                readout_mode = header["READOUT"]
            except KeyError:
                readout_mode = ""
        if readout_mode not in check_readout.split(",") or check_readout == "":
            logger.debug(
                f"Readout mode {readout_mode} not in {check_readout}. Skipping {ftsfile}."
            )

        # Binning
        try:
            x_binning = header["XBINNING"]
            y_binning = header["YBINNING"]
        except KeyError:
            x_binning = ""
            y_binning = ""
        if (
            str(x_binning) + "x" + str(y_binning) not in check_binning.split(",")
            or check_binning == ""
        ):
            logger.debug(
                f"Binning {x_binning}x{y_binning} not in {check_binning}. Skipping {ftsfile}."
            )

        # Exposure time
        try:
            exptime = header["EXPTIME"]
        except KeyError:
            try:
                exptime = header["EXPOSURE"]
            except KeyError:
                exptime = -1
        if check_exptime != "":
            for c_exp in check_exptime.split(","):
                if exptime < float(c_exp) * 0.99 or exptime > float(c_exp) * 1.01:
                    logger.debug(
                        f"Exposure time {exptime} not in {check_exptime}. Skipping {ftsfile}."
                    )
        

        # Target
        try:
            target_name = header["OBJECT"]
        except KeyError:
            try:
                target_name = header["OBJNAME"]
            except KeyError:
                try:
                    target_name = header["SOURCE"]
                except KeyError:
                    try:
                        target_name = header["TARGNAME"]
                    except KeyError:
                        target_name = ""
        if target_name not in check_target.split(",") or check_target == "":
            logger.debug(
                f"Target {target_name} not in {check_target}. Skipping {ftsfile}."
            )

        # Actual coordinates
        try:
            ra = header["OBJCTRA"]
            dec = header["OBJCTDEC"]
        except KeyError:
            try:
                ra = header["OBJRA"]
                dec = header["OBJDEC"]
            except KeyError:
                try:
                    logger.warning("No coordinates found. Using telescope coordinates.")
                    ra = header["TELRA"]
                    dec = header["TELDEC"]
                except KeyError:
                    ra = ""
                    dec = ""
        if "" not in (ra, dec):
            obj = coord.SkyCoord(ra, dec, unit="hourangle, deg")

        # Scheduled coordinates
        try:
            sched_ra = header["SCHEDRA"]
            sched_dec = header["SCHEDDEC"]
        except KeyError:
            sched_ra = ""
            sched_dec = ""
        if "" not in (sched_ra, sched_dec):
            sched_obj = coord.SkyCoord(sched_ra, sched_dec)

        if "" not in (ra, dec, sched_ra, sched_dec):
            dra = (obj.ra - sched_obj.ra).to("arcsec")
            ddec = (obj.dec - sched_obj.dec).to("arcsec")
        else:
            dra = np.nan
            ddec = np.nan

        # ZP
        try:
            zp = header["ZMAG"]
            zp_err = header["ZMAGERR"]
        except KeyError:
            zp = np.nan
            zp_err = np.nan

        # FWHM
        try:
            fwhmh = header["FWHMH"]
            fwhmhs = header["FWHMHS"]
            fwhmv = header["FWHMV"]
            fwhmvs = header["FWHMVS"]
        except KeyError:
            fwhmh = ""
            fwhmhs = ""
            fwhmv = ""
            fwhmvs = ""

        # Moon
        try:
            moon_angle = header["MOONANGL"]
            moon_phs = header["MOONPHAS"]
        except KeyError:
            moon_angle = ""
            moon_phs = ""

        print_rows.append(
            [
                ftsfile,
                target_name,
                date.jd,
                date.iso,
                filt,
                readout_mode,
                str(x_binning) + "x" + str(y_binning),
                f"{round_exptime(exptime)}",
                obj.ra.to_string(unit="hour", sep=":"),
                obj.dec.to_string(unit="deg", sep=":"),
                f"{zp:.3f}+/-{zp_err:.3f}",
                f"{fwhmh:.2f}+/-{fwhmhs:.2f}",
                f"{fwhmv:.2f}+/-{fwhmvs:.2f}",
                f"{moon_angle:.1f}",
                f"{moon_phs:.1f}",
                f"{dra:.2f}",
                f"{ddec:.2f}",
            ]
        )

    logger.debug(f"Found {len(print_rows)} files matching criteria.")

    logger.debug("Sorting by JD...")
    print_rows.sort(key=lambda x: x[2])

    table = prettytable.PrettyTable()
    table.add_rows(print_rows)
    table.set_style(prettytable.SINGLE_BORDER)
    table.field_names = [
        "FITS file",
        "Target",
        "JD",
        "UT",
        "Filter",
        "Readout",
        "Binning",
        "Exp. time [s]",
        "RA",
        "Dec",
        "ZP",
        "FWHM H [pix]",
        "FWHM V [pix]",
        "Moon angle [deg]",
        "Moon phase [%]",
        "dRA [arcsec]",
        "dDec [arcsec]",
    ]

    table.align["FITS file"] = "l"
    click.echo(table)
    click.echo()
    click.echo(f"Number of images = {len(print_rows)}")

    # fwhmh = np.median([row[11] for row in print_rows])
    # fwhmhs = np.std([row[11] for row in print_rows])
    # fwhmv = np.median([row[12] for row in print_rows])
    # fwhmvs = np.std([row[12] for row in print_rows])
    # click.echo(f"Median FWHM H = {fwhmh:.2f} +/- {fwhmhs:.2f} pix")
    # click.echo(f"Median FWHM V = {fwhmv:.2f} +/- {fwhmvs:.2f} pix")

    moon_phs = np.mean([float(row[14]) for row in print_rows])
    click.echo(f"Mean Moon phase = {moon_phs:.2f}")

    click.echo("\nMedian zero-point magnitudes")
    table = prettytable.PrettyTable()
    table.set_style(prettytable.SINGLE_BORDER)
    table.field_names = ["Readout", "Filter", "ZP"]
    for filt in ["u", "b", "g", "r", "i", "z"]:
        rows = [row for row in print_rows if row[4] == filt]
        if len(rows) == 0:
            continue
        for readout in np.unique([row[5] for row in rows]):
            zp = np.median([row[10] for row in rows if row[5] == readout])
            zp_err = np.std([row[10] for row in rows if row[5] == readout])
            table.add_row([readout, filt, f"{zp:.3f}+/-{zp_err:.3f}"])
    click.echo(table)

    return print_rows


fitslist = fitslist_cli.callback

if __name__ == '__main__':
    fitslist_cli()

