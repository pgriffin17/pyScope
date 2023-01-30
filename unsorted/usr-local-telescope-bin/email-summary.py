#!/usr/bin/env python

# Looks in image directory for recent (<24 hr old) images, creates summary email with night's statistics, emails to  observers, others
# 25 Nov 2015 RLM

# 08 Dec 2015 recognize SSON images
# 27 Jan 2016 - add pscale for new camera; skip grism images when calculating FWHM; fix FWHM calculation; use current day number to parse image list
# 14 Feb 2016 account for binning in calculating FWHM
# 18 Jan 2017 change pixel scale for CG42
# 11 OCt 2017 change pixel scale for SBIG 6303e
# 08 Feb 2018 change pixel scale for IKON L936 
# 16 Mar 2018 fix problem with UT/FWHM arrays
# 07 Apr 2019 add extinction plot, color code both plots
# 25 May 2019 change airmass correction to seeing to -0.43 based on Irwin 1966 AJ
# 10 Jun 2019 change nominal G-filter ZP (camera gain change)
# 13 Apr 2020 RLM fixed plots
# 15 Jun 2020 adjust zero-pointg-mag to 21.0 (SBIG)
# 19 May 2021 Change manager email to Caroline Roberts; Change zero-point magnitude to 21.4 (Sloan g) after cleaning corrector lens (!)
# 12 Dc 2021 Change ZPmag for AC4040 camera (RLM)

# import needed libraries
import sys, os, shutil, glob, smtplib, datetime,fnmatch
import dominate
from dominate.tags import *
import astropy.io.fits as pyfits
import datetime as dt
import numpy as np
import matplotlib as mpl
mpl.use('Agg') # Needed to use matplotlib in cron jobs
import matplotlib.pyplot as plt
plt.rcParams['axes.grid'] = True

from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart


# Define paths for images, html output file, text message
image_dir = '/usr/local/telescope/user/images/'
email_dir = '/home/talon/summary_email'                      # Folder contianing email summary daemon  and data files
html_name = '%s/observing_summary' % email_dir     # Summary date (one created each day with -daynr suffix
email_txt ='%s/email-msg.txt' % email_dir                    # General information text
email_all = '%s/gemini-emails.dat' % email_dir             # Email list  of people who always receive daily summary 

def writeline(s,face,color,n):
    hfile.write( '<font face="%s", font color ="%s", font size = %ipx>%s</font><p>' % (face,color,n,s) )
    return

def SendMail(Images, html_file, Mail_list):
    msg = MIMEMultipart()
    msg['Subject'] = 'Iowa Robotic Telescope summary: %s' % date_now()
    msg['From'] = 'robert-mutel@uiowa.edu'
    msg['To'] = 'robert-mutel@uiowa.edu'

    # Read summary file
    email_msg = open(html_file,'r')
    message = email_msg.read()
    email_msg.close()

    text = MIMEText(message,'html')
    msg.attach(text)
    if Images[0] != '':
        for image in Images:
            img_data = open(image, 'rb').read()
            image = MIMEImage(img_data, name=os.path.basename(image))
            msg.attach(image)

    s = smtplib.SMTP('localhost')
    From = 'robert-mutel@uiowa.edu'
    
    s.sendmail(From, Mail_list, msg.as_string())
    s.quit()

def date_now():
    # Construct current date string
    now = datetime.datetime.now()
    yr_now = now.year; month_now = now.month; day_now = now.day
    #str_now = '%02d%02d' % (month_now,day_now)
    #date = '%4d-%02d-%02d' % (yr_now, month_now,day_now)
    date = '%2d %s% 4d' % (day_now,now.strftime("%B"),yr_now)
    return date


# Plot FWHM vs UT
def plot_stats(UT_fwhm, FWHM, UT_zmag, D_zmag, daynr):
    
    f, ax = plt.subplots(2, sharex=True)
    f.suptitle('Zenith seeing, transparency vs. UT, %s (day %d)' % (date_now(), daynr))
    
    # Upper plot: FWHM
    colormap = plt.cm.brg_r
    normalize = mpl.colors.Normalize(vmin=2.0,vmax=5)
    ax[0].scatter(UT_fwhm, FWHM, marker = 'o',cmap=colormap,norm=normalize,  c=FWHM )
    ax[0].set_xlim(1,13)
    ax[0].set_ylim(1.5,4)
    ax[0].set_xlabel('UT hours') 
    ax[0].set_ylabel('Zenith seeing (arcsec)') 

    #Lower plot: extinction
    colormap = plt.cm.brg_r
    normalize = mpl.colors.Normalize(vmin=0,vmax=3)
    ax[1].scatter(UT_zmag, D_zmag, marker = 'o',cmap=colormap,norm=normalize, s=30,c=D_zmag)
    ax[1].set_ylim(2.0,-0.2)
    ax[1].set_ylabel('g-magnitude extinction')
    
    plt.axhline(y=0,linestyle='dashed')
    
    pltname = '%s/fwhm-ut-%s.png' % (email_dir,daynr)
    plt.savefig(pltname)
    print('Created %s' % pltname)
    return pltname


#### MAIN #####
# Plate scale (arcsec/pixel, unbinned)
pscale = 0.54 # SBIG  4040

# Nominal clear sky g-filter zero-point magnitude using  SBIG4040 
zmag_g0 = 22.2

# SSON code
sson_prefix = '7_'

# Get current day number
day_nr = dt.datetime.now().timetuple().tm_yday
day_nr =str( '%03i' % day_nr)
daynr = int(day_nr)
#print(day_nr)

# Read message file
fmsg = open(email_txt,'r')
email_msg = fmsg.read()
fmsg.close()

# Read  email_all  for email addresses of administrative receipients
f = open(email_all,'r')
email_list = [x.rstrip() for x in f.readlines()]
f.close()

# Generate list of FITS files 
fts_names = glob.glob('%s%s' % (image_dir,'*.fts'))
Obs_codes = []; Observer = []; Exp_sec = []; Date = []; UT = []; Filter = []; Source = [];FWHM = []; UT_fwhm = []
Moonphase = []; Airmass =[] ; UT_zmag = []; D_zmag = []
for fname in fts_names:
    filename = os.path.basename(fname)
    day = filename[3:6]; obscode = filename[0:3]
    if day == day_nr:
        daynr =int(day)
    else:
        continue
    # Get info from FITS header
    hdr = pyfits.getheader(fname)
    fil = hdr['FILTER'][:1]
    Filter.append(fil)
    source = hdr['OBJECT']
    Moonphase.append( hdr['MOONPHAS'])
    # Grism filter?
    if fil == '8' or fil == '9':
        grism = True
    else:
        grism = False
    
    # Clean up source name if needed (no /'s, spaces)
    source = source.replace(' ','') ; source = source.replace('/','_')
    date,ut = hdr['DATE-OBS'].split('T'); ut = ut[:-3]
    date = date.replace('-','_'); ut = ut.replace(':','')
    Date.append(date); UT.append(ut)
    Exp_sec.append(float(hdr['EXPTIME']))
    Source.append(source)
    
    # Get zmag stats from foc images
    if obscode == 'foc' and fil == 'G' and 'ZMAG' in hdr:
        ut_zmag = np.modf(hdr['JD']-0.5)[0]*24
        d_zmag = zmag_g0 - hdr['ZMAG']
        UT_zmag.append(ut_zmag); D_zmag.append(d_zmag)
    
    if fname.find(sson_prefix) > -1:
        Obs_codes.append('sso')
        Observer.append('SSON observers')
    else:
        Obs_codes.append(filename[:3])
        Observer.append(hdr['OBSERVER'])
    # Add FWHM (but skip Grism images)
    if 'FWHMH ' in hdr and not grism:
        fh = float(hdr['FWHMH ']) ; fv = float(hdr['FWHMV '])
        nbin = int(hdr['XBINNING'])
        fwhm = np.sqrt(fh*fv)*pscale * nbin
        ut_fwhm = np.modf(hdr['JD']-0.5)[0]*24
        if fwhm > 1.4:    # Don't use crazy low fwhm values
            FWHM.append(fwhm); Airmass.append(hdr['AIRMASS']); UT_fwhm.append(ut_fwhm)

# Bail if no images!
if len(Obs_codes) == 0: sys.exit('No images for %s, exiting' % date_now() )

# Calculate zenith seeing, create arrays for plotting
Z = np.array(Airmass); FZ = FWHM * (Z**-0.43)
FZ_median  =  np.median( FZ )
FZ_sigma =  np.std( FZ )

# Generate a list of unique observer codes and combined observer+administrative emails
unique_codes = list(set(Obs_codes))
unique_emails = list(set(Observer))
all_email = list(set(unique_emails + email_list))

lines = []
html_file = '%s-day%03d.html' % (html_name,daynr)
hfile = open(html_file,'w')

# Print some summary data
s = 'Iowa Robotic Telescope daily report  %s (day %i)' % (date_now(), daynr) ; writeline(s,'arial','blue',4)
hfile.write('<font face="arial">%s</font><p>' % email_msg)
s = '<b>Observing statistics</b>'; writeline(s,'arial','blue',2)
s = 'Total number of observing projects: %i' % len(unique_codes); writeline(s,'arial','black',2)
s = 'Total number of images: %i' % len(Obs_codes); writeline(s,'arial','black',2)
s =  'Moon phase range: %.0f%% - %.0f%%' % ( min(Moonphase), max(Moonphase) ); writeline(s,'arial','black',2)
s = 'Median zenith seeing: %.1f +/- %.1f arcsec'  % (FZ_median,FZ_sigma); writeline(s,'arial','black',2)
Exp = []; sources = []
for code in unique_codes:
    ncode = Obs_codes.count(code)
    indices = [i for i, x in enumerate(Obs_codes) if x == code]
    sources.append(list ( set ( [Source[i] for i in indices] )))
    Exp_sum = int(sum([Exp_sec[i] for i in indices]) /60.)
    obs = Observer[indices[0]].split('@')[0]
    ut = np.array( [int(float(x)) for x in  [UT[i] for i in indices] ])
    utmin = np.min(ut); utmax = np.max(ut)  
    Exp.append(Exp_sum)
    line = [code, obs, ncode,'%i' % Exp_sum, '%04d - %04d' % (utmin,utmax) ]
    lines.append(line)

# Make a table of statistics, one line per observer code
s = 'Total observing time: %.1f hrs' % (sum(Exp)/60.); writeline(s,'arial','black',2)
hfile.write('<font face="arial">')
h = html()
col_align=['center','left','right', 'right','center']
col_style=['','','',' ',' ']
with h.add(table(cellpadding = "4", style = "border: 1px solid #000000; border-collapse: collapse;", border = "1")):
    r = tr()
    [r.add(th(x)) for x in ['Code', 'Observer', 'N image', 'Time (min)', 'UT range']]
    for line in lines:
        r = tr()
        [r.add(td(line[i], align = col_align[i], style = col_style[i])) for i in range(5)]

htmlcode = str(h)[9:-8]
hfile.write(htmlcode + '</font><p>\n')

# list all observed objects by observer code
s = '<b>Objects observed by observer code</b>'; writeline(s,'arial','blue',2)
for n in range(len(unique_codes)):
    code = unique_codes[n]; src_list = ', '.join(sources[n])
    s =  '<b>%s</b>: %s' % (code, src_list) ; writeline(s,'arial','black',2)

# Supply email for list removal
s = 'You are receiving this email either because you had completed observations or you are on the telescope notification list. If you wish to be removed from this list, or have any questions or comments, please email: '
writeline(s,'arial','black',2)
h = html()
with h:
    p('Caroline Roberts AT uiowa.edu', href = 'mailto:caroline-roberts@uiowa.edu')
htmlcode = str(h)[9:-8]
#print(htmlcode)
#htmlcode=HTML.link('Robert-Mutel AT uiowa.edu', 'mailto:robert-mutel@uiowa.edu')
hfile.write(htmlcode)
hfile.close()

#  FWHM vs UT plot
pltname = plot_stats(UT_fwhm, FZ,  UT_zmag, D_zmag, daynr)

# Send mail, attaching images if desired (stub for now)
Images = [pltname]
#Mail_list = ['jua.kamp451@gmail.com']
#Mail_list = ['rmutel@gmail.com']
Mail_list = all_email
SendMail(Images, html_file, Mail_list)
print('Emailed %s to users: %s' % (html_file,Mail_list))
