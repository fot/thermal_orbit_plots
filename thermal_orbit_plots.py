import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
from copy import deepcopy
import pprint
from os.path import expanduser

import Ska
import Ska.engarchive.fetch_eng as fetch
from Chandra.Time import DateTime

home = expanduser("~")
sys.path.append(home + '/AXAFLIB/pylimmon/')

# import gretafun
import pylimmon

plt.rcParams['xtick.major.pad'] = '10'


def importcolors():
    ''' Return a dictionary of colors along with their hex and RGB values'''

    filename = './colors.csv'
    fields = ['color','hex','r','g','b']
    colorarray = np.genfromtxt(filename, delimiter=',', dtype=None, names=fields,
                           comments=None)
    colordict = {}
    
    for (name,hex,r,g,b) in colorarray:
        
        if r > 1 or g > 1 or b > 1:
            r = r/255.0
            g = g/255.0
            b = b/255.0
            
        colordict[name] = {'hex':hex, 'rgb':[r,g,b]}

    return colordict

colordict = importcolors()


#-------------------------------------------------------------------------------------------------
# Originally from greta_parse.py
#-------------------------------------------------------------------------------------------------


def parsedecplot(decfile, removewidechars=True):
    '''Parse a GRETA dec plot file to extract plotting data. This will not
    grab text display data.
    '''

    def finddecstring(pattern, string, rtype='string', split=False):
        '''Search a single line in a GRETA dec file for plotting
        information.

        "pattern" is the name of the GRETA keyword
        "string" is the string to be searched
        "rtype" is the returned data type (currently string, int, or float)
        "split" is a flag to allow the user to request the returned data be
            split into a list
        '''

        rtype = str(rtype).lower()
        p1 = re.compile('^' + pattern + '\s+(.*)$', re.MULTILINE)
        r = p1.search(string)
        if not isinstance(r, type(None)):
            rval = r.group(1)
            if split:
                rval = rval.split()
            if rtype == 'string':
                if split:
                    rval = [val.strip() for val in rval]
                else:
                    rval = rval.strip()
            if rtype=='float':
                if split:
                    rval = [float(n) for n in rval]
                else:
                    rval = float(rval)
            if rtype=='int':
                if split:
                    rval = [int(n) for n in rval]
                else:
                    rval = int(rval)
            return rval

    #filename = 'Orbit_Plots/T_STT72_ISIM_ACIS.dec'
    infile = open(decfile,'rb')
    body = infile.read()
    infile.close()

    decplots = {}

    decplots['DTITLE'] = finddecstring('DTITLE', body)
    decplots['DSUBTITLE'] = finddecstring('DSUBTITLE', body)
    decplots['DTYPE'] = finddecstring('DTYPE', body, split=True)
    decplots['DTYPE'][1] = int(decplots['DTYPE'][1])
    xdata = finddecstring('DXAXIS', body, rtype='float', split=True)
    decplots['DXAXIS'] = [60*d for d in xdata]

    r = re.split('\nPINDEX',body) # Newline excludes commented out plots
    decplots['numplots'] = len(r)-1
    plots = {}
    for plotdef in r[1:]:
        num = int(re.match('\s+(\d+)[.\n]*', plotdef).group(1))
        plots[num] = {}
        plots[num]['PINDEX'] = num
        plots[num]['PTRACES'] = finddecstring('PTRACES', plotdef,
                                              rtype='int')
        plots[num]['PBILEVELS'] = finddecstring('PBILEVELS', plotdef,
                                                rtype='int')

        plots[num]['PTITLE'] = finddecstring('PTITLE', plotdef)
        plots[num]['PYLABEL'] = finddecstring('PYLABEL', plotdef)
        plots[num]['PGRID'] = finddecstring('PGRID', plotdef, rtype='int')
        plots[num]['PLEGEND'] = finddecstring('PLEGEND', plotdef,
                                              rtype='int')
        plots[num]['PYAXIS'] = finddecstring('PYAXIS', plotdef,
                                             rtype='float', split=True)
        plots[num]['PYAUTO'] = finddecstring('PYAUTO', plotdef, rtype='int')

        # cycle through all TINDEX in a similar way to how you cycle through
        # PINDEX
        t = re.split('\nTINDEX',plotdef) # Newline excludes commented out msids
        traces = {}
        for tracedef in t[1:]:
            tnum = int(re.match('\s+(\d+)[.\n]*', tracedef).group(1))
            traces[tnum] = {}
            traces[tnum]['TINDEX'] = tnum
            traces[tnum]['TMSID'] = finddecstring('TMSID', tracedef)
            traces[tnum]['TNAME'] = finddecstring('TNAME', tracedef)
            traces[tnum]['TCOLOR'] = finddecstring('TCOLOR', tracedef)
            traces[tnum]['TCALC'] = finddecstring('TCALC', tracedef)
            traces[tnum]['TSTAT'] = finddecstring('TSTAT', tracedef)
            if removewidechars:
                if traces[tnum]['TMSID'] != None:
                    traces[tnum]['TMSID'] = traces[tnum]['TMSID'].replace('_WIDE', '')
                    traces[tnum]['TMSID'] = traces[tnum]['TMSID'].replace('_wide', '')
        plots[num]['traces'] = traces

        # TBLINDEX - do this separate in case they are intermingled with
        # TINDEX definitions
        tb = re.split('\nTBLINDEX',tracedef) # Newline excludes comments
        if len(tb) > 1:
            tbtraces = {}
            for tbtracedef in tb[1:]:
                tbnum = re.match('\s+(\d+)[.\n]*', tbtracedef)
                tbnum = int(tbnum.group(1))
                tbtraces[tbnum] = {}
                tbtraces[tbnum]['TBINDEX'] = tbnum
                tbtraces[tbnum]['TMSID'] = finddecstring('TMSID', tbtracedef)
                tbtraces[tbnum]['TNAME'] = finddecstring('TNAME', tbtracedef)
                tbtraces[tbnum]['TCOLOR'] = finddecstring('TCOLOR', tbtracedef)
                if removewidechars:
                    if traces[tnum]['TMSID'] != None:
                        traces[tnum]['TMSID'] = traces[tnum]['TMSID'].replace('_WIDE', '')
                        traces[tnum]['TMSID'] = traces[tnum]['TMSID'].replace('_wide', '')
            plots[num]['tbtraces'] = tbtraces

    decplots['plots'] = plots

    return decplots


#-------------------------------------------------------------------------------------------------
# Original dechelper.py code
#-------------------------------------------------------------------------------------------------

class fetchobject(Ska.engarchive.fetch.Msid):
    def __init__(self, msid, times, vals, tstart, tstop, stat=None):
        self.times = times # This gets overwritten if stats are requested
        self.vals = vals
        self.MSID = 'DP_' + msid.upper()
        self.datestart = DateTime(tstart).date
        self.datestop = DateTime(tstop).date
        self.tstart = DateTime(tstart).secs
        self.tstop = DateTime(tstop).secs
        if stat:
            #import code
            #code.interact(local=locals())
            stats = getstats(self, stat)
            self.__dict__.update(stats)


def getstats(data, stat):

    tstart = data.times[0]
    tstop = data.times[-1]
    dt = np.diff(data.times)[0]

    if stat.lower() == '5min':
        pts = int(round(328 / dt))
    elif stat.lower() == 'daily':
        pts = int(round(86400 / dt))

    if stat:
        # Note this assumes no missing values...
        numintervals = int(len(data.times) / pts)

        timemat = np.reshape(data.times[:numintervals * pts], (numintervals, pts))
        valmat = np.reshape(data.vals[:numintervals * pts], (numintervals, pts))

        meantimes = np.mean(timemat, axis=1)
        meanvals = np.mean(valmat, axis=1)
        minvals = np.min(valmat, axis=1)       
        maxvals = np.max(valmat, axis=1)
        stds = np.std(valmat, axis=1)
        midvals = valmat[:, round(pts/2)]

        returndict = {'times':meantimes, 'means':meanvals, 'mins':minvals,
                      'maxes':maxvals, 'midvals':midvals, 'stds':stds}
        return returndict

def alt_msid_fetch(msid, time1, time2, stat=None):

    if 'OBAHCHK' in msid.upper():

        msids = ['OOBTHR08', 'OOBTHR09', 'OOBTHR10', 'OOBTHR11', 'OOBTHR12', 
                 'OOBTHR13', 'OOBTHR14', 'OOBTHR15', 'OOBTHR17', 'OOBTHR18', 
                 'OOBTHR19', 'OOBTHR20', 'OOBTHR21', 'OOBTHR22', 'OOBTHR23', 
                 'OOBTHR24', 'OOBTHR25', 'OOBTHR26', 'OOBTHR27', 'OOBTHR28', 
                 'OOBTHR29', 'OOBTHR30', 'OOBTHR31', 'OOBTHR33', 'OOBTHR34', 
                 'OOBTHR35', 'OOBTHR36', 'OOBTHR37', 'OOBTHR38', 'OOBTHR39', 
                 'OOBTHR40', 'OOBTHR41', 'OOBTHR42', 'OOBTHR45', 'OOBTHR46']

        data = fetch.Msidset(msids, time1, time2)
        data.interpolate()

        maxes = data[msids[0]].vals
        mins = data[msids[0]].vals
        for name in data.keys():
            maxes = np.max((maxes, data[name].vals), axis=0)
            mins = np.min((mins, data[name].vals), axis=0)
            

        return fetchobject('OBAHCHK', data.times, maxes-mins, time1, time2, stat)


    elif 'HADG' in msid.upper():

        msids = ['OHRMGRD3', 'OHRMGRD6']

        data = fetch.Msidset(msids, time1, time2)
        data.interpolate()

        vals = np.max((data[msids[0]].vals, data[msids[1]].vals), axis=0)

        return fetchobject('HADG', data.times, vals, time1, time2, stat)


    elif 'POBA' in msid.upper():

        data = fetch.Msid('DP_POBAT', time1, time2)

        return fetchobject('POBA', data.times, data.vals, time1, time2, stat)


    elif 'PSUM' in msid.upper():

        data = fetch.Msid('DP_PABH', time1, time2)

        return fetchobject('PSUM', data.times, data.vals, time1, time2, stat)


    elif 'DUTYCYCLE' in msid.upper():

        msids = ['4OHTRZ53', '4OHTRZ54', '4OHTRZ55', '4OHTRZ57']

        data = fetch.Msidset(msids, time1, time2)
        data.interpolate()

        RTOTAL = 1./( 1/94.1 + 1/124.3 + 1/126.8 + 1/142.3)
        DC1 = np.abs(data['4OHTRZ53'].raw_vals)/94.1 + np.abs(data['4OHTRZ54'].raw_vals)/124.3
        DC2 = np.abs(data['4OHTRZ55'].raw_vals)/126.8 + np.abs(data['4OHTRZ57'].raw_vals)/142.3
        DUTYCYCLE = RTOTAL*(DC1+DC2)

        return fetchobject('DUTYCYCLE', data.times, DUTYCYCLE, time1, time2, stat)




#-------------------------------------------------------------------------------------------------
# Original decplotdata.py code
#-------------------------------------------------------------------------------------------------

class plotdata(object):
    ''' Retrieve, configure, and return requested telemetry.

    This class/set of routines are required for these reasons:

    1) The plotting function this was built for requested a particular stat 
       (based on given data), which may not be available in the archive.

    2) State based telemetry does not include 5min or daily stats (bilevels).

    3) The plotting function expects to see "plotdata" and "plottimes" 
       attributes.

    4) The plotting functions expects to see limit information if available.

    5) Some "MSIDs" do not exist in the archive, but can be calculated. In
       this case, a simulated fetch.Msid class is returned.

    6) There are some time periods that should be removed from the data to 
       make the plots more relevant to monitoring trends during nominal
       operations.

    '''

    def __init__(self, msid, tstart, tstop, plotstat=None, fetchstat=None, type='numeric'):
        self.msid = msid
        self.time1 = tstart
        self.time2 = tstop
        self.plotstat = plotstat # calculated stat as defined in dec file
        self.fetchstat = fetchstat # archive stat ('5min', 'daily', or None)
        self.__dict__.update(self._gettracetelemetry())
        self.__dict__.update(self._getploplotstats())
        
    def _gettracetelemetry(self):
        '''Retrieve, configure and return telemetry.
        '''

        telem = self._fetchhelper(stat=self.fetchstat)

        # This ensures plotstat is defined as a string to make the tests below
        # more straightforward.
        if isinstance(self.plotstat,type(None)):
            self.plotstat = ''

        # Define name to use, this is used so labels can be added without
        # modifying the original msid name
        name = self.msid

        # Grab/calculate the requested data
        if self.plotstat.lower() == 'max':
            # It is assumed that all msids that use this stat are numeric
            if self.fetchstat:
                telem.plotdata = telem.maxes
            else:
                telem.plotdata = np.maximum.accumulate(telem.vals)

            telem.plottimes = telem.times

            name = name + '-max'


        elif self.plotstat.lower() == 'min':
            # It is assumed that all msids that use this stat are numeric
            if self.fetchstat:
                telem.plotdata = telem.mins
            else:
                telem.plotdata = np.minimum.accumulate(telem.vals)

            telem.plottimes = telem.times

            name = name + '-min'


        elif self.plotstat.lower() == 'mean':

            # Just look to see if mean values are present, if not then
            # calculate daily values from the full resolution data.

            if 'means' in telem.__dict__.keys():
                # Only numeric telemetry statistics objects will include means.
                telem.plotdata = telem.means
                telem.plottimes = telem.times

            else:
                # Since heater on/off data is stored as a state value,
                # statistics for these data are not calculated and stored
                # in the engineering archive.
                #
                # Remember that the character values are replaced with their 
                # corresponding raw values in the self._gettelemetry function.
                #
                # Since state based data will require the stats to be 
                # recalculated from the original data, refetch the full 
                # resolution data here and then calculate the means.
                telem = self._fetchhelper(stat=None)
                datadict = self._gettelemstats(telem, stat='daily')

                # Here you are copying over all the new stats, but you only
                # need the means and times
                telem.__dict__.update(datadict)
                telem.plotdata = telem.means
                telem.plottimes = telem.times

            name = name + '-avg'

            
        else:
            telem.plotdata = telem.vals
            telem.plottimes = telem.times

        try:
            safetylimits = pylimmon.get_safety_limits(name)
            # safetylimits = gretafun.getSafetyLimits(telem)
        except:
            safetylimits = {}

        try:
            glimmonlimits = pylimmon.get_latest_glimmon_limits(name)
            # glimmonlimits = gretafun.getGLIMMONLimits(name)
        except KeyError:
            glimmonlimits = {}


        tracedata = {'msid': self.msid,
                     'name': name,
                     'telem':telem,
                     'safetylimits':safetylimits,
                     'glimmonlimits':glimmonlimits,
                     }

        return tracedata


    def _fetchhelper(self, stat):
        '''Fetch telemetry.

        This implements a set of helper functions that can be used to
        calculate or fetch from greta data that are not in the engineering
        archive. 
        '''

        msid = self.msid.strip() # Sometimes there is a trailing space
        time1 = DateTime(self.time1).secs - 5 * 24 * 3600
        time2 = DateTime(self.time2).secs + 5 * 24 * 3600

        # Much of the code that handles the telemetry is built to expect these
        # fields, even if they are empty. Otherwise there would have to be
        # tests to see if telem exists or not everywhere.
        class emptyfetchobject(Ska.engarchive.fetch.Msid):
            def __init__(self):
                self.times = np.array([])
                self.vals = np.array([])
                self.MSID = ''
                self.datestart = ''
                self.datestop = ''
                self.tstart = ''
                self.tstop = ''

        telem = emptyfetchobject()
       
        try:
            # Try to fetch from the engineering archive
            telem = fetch.Msid(msid, time1, time2, stat=stat)
            
        except ValueError, e:

            print e

            try:
                if stat:
                    stat = "'" + str(stat) + "'"
                else:
                    stat = None

                print('Trying alternate MSID definition for: {}'.format(msid.upper()))
                telem = alt_msid_fetch(msid.upper(), time1, time2, statname=stat)

            except:
                # Then this msid is probably not implemented, move on with
                # the other plots
                #debug(locals())
                print('%s not in engineering archive, and not defined ' \
                      'elsewhere.'%msid)

        
        if any(telem.times):

            # Remove unwanted data, if data exists
            telem = self._removetimescaller(telem, stat)                   

            # Replace character values with their raw values for state based msids
            telem.type = 'numeric'

            if isinstance(telem.vals[0], type('')):
                telem.vals = np.abs(telem.raw_vals)
                telem.type = 'state'

        return telem


    def _removetimescaller(self, telem, stat):
        # There should be a better way to remove bad/unwanted data, but this
        # works for now

        # First remove the padded data
        telem = self._removetimes(telem, stat, 
                                DateTime(self.time1).secs - 10 * 24 * 3600,
                                DateTime(self.time1).secs)
        telem = self._removetimes(telem, stat, DateTime(self.time2).secs,
                                DateTime(self.time2).secs + 10 * 24 * 3600)
        
        telem = self._removetimes(telem, stat, '2011:149:00:00:00',
                                  '2011:153:00:00:00')
        telem = self._removetimes(telem, stat, '2011:186:00:00:00',
                                  '2011:195:00:00:00')
        telem = self._removetimes(telem, stat, '2011:299:00:00:00',
                                  '2011:306:00:00:00')
        telem = self._removetimes(telem, stat, '2012:149:00:00:00',
                                  '2012:152:00:00:00')
        return telem


    def _removetimes(self, telem, stat, t1, t2):
        
        # Remove 2011 Safe Mode Data and October NSM
        ind1 = telem.times < DateTime(t1).secs
        ind2 = telem.times > DateTime(t2).secs
        ind = ind1 | ind2
     
        telem.times = telem.times[ind]
        telem.vals = telem.vals[ind]
        if stat:
            #if hasattr(telem, 'maxes'):
            if 'maxes' in telem.__dict__.keys():
                telem.maxes = telem.maxes[ind]
                telem.mins = telem.mins[ind]
                telem.means = telem.means[ind]
                telem.midvals = telem.midvals[ind]

        return telem

    def _gettelemstats(self, data, stat):

        tstart = data.times[0]
        tstop = data.times[-1]
        dt = np.mean(np.double(np.diff(data.times)))

        if stat.lower() == '5min':
            pts = int(round(328 / dt))
        elif stat.lower() == 'daily':
            pts = int(round(86400 / dt))

        if stat:
            # Note this assumes no missing values...
            numintervals = int(len(data.times) / pts)

            timemat = np.reshape(data.times[:numintervals * pts], 
                                 (numintervals, pts))
            valmat = np.reshape(data.vals[:numintervals * pts], 
                                (numintervals, pts))

            meantimes = np.mean(timemat, axis=1)
            meanvals = np.mean(valmat, axis=1)
            minvals = np.min(valmat, axis=1)       
            maxvals = np.max(valmat, axis=1)
            stds = np.std(valmat, axis=1)
            midvals = valmat[:, round(pts/2)]

            returndict = {'times':meantimes, 'means':meanvals, 'mins':minvals,
                          'maxes':maxvals, 'midvals':midvals, 'stds':stds}
            return returndict


    def _getploplotstats(self):
        if 'maxes' in self.telem.__dict__.keys():
            maxval = np.max(self.telem.maxes)
            minval = np.min(self.telem.mins)

        else:
            maxval = np.max(self.telem.vals)
            minval = np.min(self.telem.vals)

        plotmax = np.max(self.telem.plotdata)
        plotmin = np.min(self.telem.plotdata)
        plotmean = np.mean(np.double(self.telem.plotdata))
        plotstd = np.std(np.double(self.telem.plotdata))

        meanval = np.mean(np.double(self.telem.vals))
        stdval = np.std(np.double(self.telem.vals))
        mintime = np.min(self.telem.times)
        maxtime = np.max(self.telem.times)

        # if there is TDB info, throw this in as well
        if 'tdb' in self.telem.__dict__.keys():
            tdb = self.telem.tdb
        else:
            tdb = None

        return {'plotmax':plotmax, 'plotmin':plotmin, 'plotmean':plotmean, 
                'telemmax':maxval, 'telemmin':minval, 'telemmean':meanval,
                'plotstd':plotstd, 'telemstd':stdval, 'mintime':mintime, 
                'maxtime':maxtime, 'tdb':tdb}





#-------------------------------------------------------------------------------------------------
# Original pydecplot.py code
#-------------------------------------------------------------------------------------------------


def debug(locals):
    import code
    code.interact(local=locals)

class plotdec(object):
    ''' Use the GRETA dec file to produce enhanced orbit plots.

    Add more description here

    Time period must be more than one day, preferably more than two, due to
    how daily stats are calculated (specifically daily means)
    '''

    def __init__(self, decfile, time1=None, time2=None, fgcolor=[1,1,1],
                 bgcolor=[0.15, 0.15, 0.15], orientation='landscape',
                 plotltt=False):

        self.decfile = decfile
        self.decplots = parsedecplot(self.decfile)
        self.colors = importcolors()
        self.plotltt = plotltt
        self.ltttimespan = 60 * 24 * 3600
        self.defaultcolors = ['RED', 'YELLOW', 'GREEN', 'AQUA', 'PINK', 'WHEAT',
                              'GREY', 'BROWN', 'BLUE', 'BLUEVIOLET', 'CYAN',
                              'TURQUIOSE', 'MAGENTA', 'SALMON', 'WHITE']
        self.plotinfo = {'bgcolor':bgcolor,
                         'fgcolor':fgcolor,
                         'top':0.88,
                         'wspace':None,
                         'hspace':0.3,
                         'binarylocation':[None]*self.decplots['numplots']}
        if orientation.lower() == 'landscape':
            sizeparam =  {'width':18,
                          'height':10,
                          'left':0.05,
                          'bottom':0.05,
                          'right':0.78,
                          'stampfontsize':10,
                          'datefontsize':10,
                          'statsfontsize':8,
                          'statsvscalefact':1.3,
                          'lttspace':0.2}
        else:
            sizeparam =  {'width':8.5,
                          'height':11,
                          'left':0.07,
                          'bottom':0.05,
                          'right':0.69,
                          'stampfontsize':6,
                          'datefontsize':6,
                          'statsfontsize':6,
                          'statsvscalefact':1.1,
                          'lttspace':0.2}
        self.plotinfo.update(sizeparam)
        self.plotinfo['location'] = self._getplotloc()
        self._addbinaryplotloc()

        if plotltt:
            self._addlttsplots()
            
        if time2:
            self.time2 = DateTime(time2).date
        else:
            self.time2 = DateTime().date

        if time1:
            self.time1 = DateTime(time1).date
        else:
            self.time1 = DateTime(DateTime(time2).secs - 10*24*3600).date

        self._plotfigure()
        plt.draw()

    def _addlttsplots(self):
        ''' Make space for LTT plots on the left.

        This fills out the location array for all LTT plots based on the 
        number of plots.
        '''
        
        self.plotinfo['lttslocation'] = deepcopy(self.plotinfo['location'])
        for loc in self.plotinfo['lttslocation']:
            loc[0] = self.plotinfo['left']
            loc[2] = self.plotinfo['lttspace']
        

    def _getplotloc(self):
        ''' Generate the location and sizing for all primary plots.

        This does not include LTT and binary plots, Binary and LTT plot sizing
        are defined in another function. This function will shift the primary 
        plots to make room for LTT plots if LTT plots are requested.
        '''
        numplots = self.decplots['numplots']
        hspace = self.plotinfo['hspace']
        top = self.plotinfo['top']
        bottom = self.plotinfo['bottom']
        wspace = self.plotinfo['wspace']
        right = self.plotinfo['right']

        if self.plotltt:
            left = self.plotinfo['left'] + self.plotinfo['lttspace'] + 0.05
        else:
            left = self.plotinfo['left']
        
        # Spacing between each plot is a fraction of the available space for
        # each plot
        spacing = hspace * (top - bottom) / numplots

        # Total available space to plot stuff / number of plots
        plotspace = (top - bottom - spacing * (numplots - 1)) / numplots

        plotloc = []
        for n in range(numplots - 1, -1, -1):
            # Origin and size of each plot
            # [x, y, width, height]
            plotloc.append([left, bottom + plotspace * n + spacing * n,
                            right - left, plotspace])

        return plotloc
    

    def _addbinaryplotloc(self):
        ''' Define location and sizing for binary plots.

        This is intended to mimic the binary plotting capabilities included 
        in GRETA. The addition of binary data should decrease the size of the 
        primary plot.
        '''

        # these two are in the figure coordinate system
        baseheight = 0.01 # Height without any traces yet
        traceheight = 0.007 # Amount to add for each trace

        # Step through each plot
        for p in self.decplots['plots'].keys():

            # if the binary plot key doesn't have a NoneType
            if self.decplots['plots'][p]['PBILEVELS']:
                
                # This is the number of binary traces for a particular plot
                # You could get this number either from PBILEVELS or just
                # counting the number of keys as is done below.
                numtbtraces = len(self.decplots['plots'][p]['tbtraces'].keys())

                # Total height for the binary trace plot (fits under regular
                # plot, which is required)
                totalheight = baseheight + traceheight * numtbtraces

                # Copy over the primary plot parameters
                binaryloc = deepcopy(self.plotinfo['location'][p])

                # Only need to change the height, and then make room by
                # modifying/shrinking the primary plot
                binaryloc[3] = totalheight
                self.plotinfo['location'][p][1] = \
                                  self.plotinfo['location'][p][1] + totalheight
                self.plotinfo['location'][p][3] = \
                                  self.plotinfo['location'][p][3] - totalheight

                # Copy over the new binary plot parameters
                self.plotinfo['binarylocation'][p] = binaryloc


    

    def _getcolor(self, tcolor=None, tracenum=0):
        if isinstance(tcolor, type(None)):
            tcolor = self.defaultcolors[tracenum]
        elif tcolor.lower() in self.colors.keys():
            tcolor = self.colors[tcolor.lower()]['rgb']
        else:
            tcolor = self.defaultcolors[tracenum]
        return tcolor


    def _generategrid(self, ax):
        #ax.xaxis.set_minor_locator(AutoMinorLocator())
        #ax.yaxis.set_minor_locator(AutoMinorLocator())
        #ax.tick_params(axis='both', which='minor', length=6, color=[0.8, 0.8, 0.8])
        ax.grid(b=True, which='both', color=self.plotinfo['fgcolor'])
        ax.set_axis_bgcolor(self.plotinfo['bgcolor'])
        ax.set_axisbelow(True)

    
    def _generatelegend(self, ax, plotnum):
        bbox_to_anchor = (0, 1, 1, 0.03)

        numtraces = len(self.decplots['plots'][plotnum]['traces'].keys())
        leg = ax.legend(bbox_to_anchor=bbox_to_anchor, loc=8, ncol=numtraces,
                                     mode="expand", prop= {'size':10},
                                     borderaxespad=0.,handletextpad=0.1)
        frame  = leg.get_frame()  
        frame.set_facecolor(self.plotinfo['bgcolor'])

        for t in leg.get_texts():
            t.set_fontsize(self.plotinfo['statsfontsize'])
            t.set_color(self.plotinfo['fgcolor'])


    def _configurexaxis(self, ax, plotnum, t1=None, t2=None, numticks=11):
        ''' Define x axis ticks.

        Only plots tick marks, date labels are plotted separately

        '''

        if not t1:
            t1 = DateTime(self.time1).secs

        if not t2:
            t2 = DateTime(self.time2).secs

        # Generate a list of tick mark locations, the corresponding labels will
        # empty strings, the actual date labels are only plotted at the bottom
        # in the figure coordinate system, not the axis coordinate system.
        xtik = np.linspace(t1, t2, numticks)

        xlab = ['']

        ax.set_xticks(xtik)
        ax.set_xticklabels(xlab)

        ax.set_xlim(t1, t2)


    def _plotdatelabels(self, plotloc, t1=None, t2=None, numticks=11, 
                        lttplot=False):
        ''' Define x axis labels/ticks.

        This plots the time stamps manually as annotations in the figure 
        coordinate system. This should be plotted in the figure coordinate
        system to avoid conflicts with binary plots. If the bottom plot
        includes a binary plot, the primary plot axis can't be used to show
        the axis labels (the binary plot will cover these labels). Rather 
        than writing the logic to place the labels in one plot function or 
        another, and to remember the right vertical location for plotting
        the LTT date labels, the date label function is kept separate.

        '''

        if not t1:
            t1 = DateTime(self.time1).secs

        if not t2:
            t2 = DateTime(self.time2).secs

        # Vertical location of date labels
        yloc = plotloc[1] - 0.01

        # Generate a list of tick mark locations, the corresponding labels will
        # either include date stamps, or empty strings so that only the bottom 
        # plot shows the dates.
        xtik = np.linspace(t1, t2, numticks)

        # Generate date labels
        xlabels = DateTime(xtik).date
        xlabels = [name[:8] + ' \n' + name[9:17] + ' ' for name in xlabels]

        # Generate version of xtik in axis coordinate system
        xtik_ax = np.linspace(0, 1, numticks)

        # Convert axes x coordinate to figure coordinate system
        #
        # new x = axes ratio / axes width in fig coord. + axes left coord.
        if lttplot:
            xtik_fig = [x_ax * self.plotinfo['lttspace'] + 
                        self.plotinfo['left'] for x_ax in xtik_ax]
        else:
            xtik_fig = [x_ax * plotloc[2] + plotloc[0] for x_ax in xtik_ax]
            

        for xloc, xlab in zip(xtik_fig, xlabels):
            self.fig.text(xloc, yloc, xlab, ha="center", va="top",
                          size=self.plotinfo['datefontsize'],
                          color=self.plotinfo['fgcolor'])


   
    def _configureyaxis(self, ax, tracestats=None):
       
        ax.tick_params(axis='y', labelsize=10,
                       labelcolor=self.plotinfo['fgcolor'],
                       color=self.plotinfo['fgcolor'])
        if tracestats:

            miny = 1e6
            maxy = -1e6
            for statkey in tracestats:
                miny = np.min((miny, tracestats[statkey].plotmin))
                maxy = np.max((maxy, tracestats[statkey].plotmax))
        else:
            miny = 0
            maxy = 1
        
        miny = miny - (maxy - miny) * 0.05
        maxy = maxy + (maxy - miny) * 0.05
        ax.set_ylim(miny, maxy)
        

    def _writestats(self, tracestats, plotnum):

        # Find longest msid name string
        longest = 9 # minimum width
        for stat in tracestats.keys():
            longest = np.max((longest, len(tracestats[stat].name[4:])))

        msidpad = longest + 1
        stringconstructor = '%' + str(msidpad) + 's%9s%9s%10s%10s\n'

        text = stringconstructor%('', 'Trending', 'Trending', '', 'Max This')
        text = text + stringconstructor%('Trace', 'Caution', 'Caution', 'Prior',
                                             'Time')
        text = text + stringconstructor%('Name', 'Low', 'High', 'Year Max',
                                             'Period')
        text = text + '\n'

        for name in tracestats.keys():
            if name[:3] == 'stt':
                if 'glimmonlimits' in tracestats[name].__dict__.keys():
                    if tracestats[name].glimmonlimits is not None:
                        if tracestats[name].glimmonlimits.has_key('caution_low'):
                            cautionlow = str(tracestats[name].glimmonlimits\
                                             ['caution_low'])
                            cautionhigh = str(tracestats[name].glimmonlimits\
                                              ['caution_high'])
                    else:
                        cautionlow = 'None'
                        cautionhigh = 'None'                    
                else:
                    cautionlow = 'None'
                    cautionhigh = 'None'

                # Most yearly data will have maxes, since only the daily statistics
                # are fetched (to save memory). Some yearly data, such as heater
                # state data need to be queried at full resolution so that the
                # daily statistics can be manually calculated. If this data is
                # manually calculated, then the data is not stored in the maxes,
                # instead it is merely stored in the yeartelem.plotdata attribute.
                lttname = 'ltt' + name[3:]
                if lttname in tracestats.keys():
                    yearmax = str('%10.5f'%tracestats[lttname].telemmax)
                else:
                    yearmax = 'None'

                # Recent data is not based on statistics so it will not likely have
                # a 'maxes' attribute
                recentmax = tracestats[name].telemmax   
                maxplotted = str('%10.5f'%recentmax)
                
                text = text + stringconstructor%(name[4:], cautionlow, cautionhigh,
                                                 yearmax, maxplotted)

        xloc = self.plotinfo['right'] + 0.01
        yloc = self.plotinfo['location'][plotnum][1] + \
               self.plotinfo['location'][plotnum][3] + 0.02
        
        stats = self.fig.text(xloc, yloc, text, ha="left", va="top",
                             size=self.plotinfo['statsfontsize'],
                             family='monospace', color=self.plotinfo['fgcolor'])


    def _getplotdata(self, plotnum, tracenum, time1=None, 
                     fetchstat=None, plotstat=None, binaryplot=False):

        # Figure out what the MSID name is
        if binaryplot:
            tmsid = self.decplots['plots'][plotnum]['tbtraces'][tracenum]['TMSID']
            tcalc = None
            tstat = None
        else:
            tmsid = self.decplots['plots'][plotnum]['traces'][tracenum]['TMSID']
            tcalc = self.decplots['plots'][plotnum]['traces'][tracenum]['TCALC']
            tstat = self.decplots['plots'][plotnum]['traces'][tracenum]['TSTAT']

        # If no initial time is specified, use the initial time passed to 
        # this class. Remember you can't use self to initialize another 
        # keyword in the function definition above.
        if not time1:
            time1 = self.time1

        if tcalc:
            name = tcalc
        elif tmsid:
            name = tmsid
        else:
            raise ValueError("MSID name missing")

        # if plotstat is currently None, and the decfile included a tstat, use
        # the stat calculation requested by tstat
        if (not plotstat) and tstat:
            plotstat=tstat

        telem = None
        tracedata = None
        try:
            tracedata = plotdata(name, time1, self.time2, 
                                             plotstat=tstat, 
                                             fetchstat=fetchstat)
            telem = tracedata.telem
            del tracedata.telem

        except ValueError as e:
            # debug(locals())
            print('Empty plot.')
            print('Returned Error (in self._plotaxis): %s'%e)

        return telem, tracedata


    def _plotaxis(self, plotnum):

        # Define attributes to keywords to simplify
        bgcolor = self.plotinfo['bgcolor']
        fgcolor = self.plotinfo['fgcolor']
        ax = self.axlist[plotnum]

        # Add title and other plot info at this level
        ax.set_title(self.decplots['plots'][plotnum]['PTITLE'], fontsize=12,
                     color=fgcolor, position=[0.5, 1.15])


        tracestats = {}
        for tracenum in self.decplots['plots'][plotnum]['traces'].keys():

            telem, tracedata = self._getplotdata(plotnum, tracenum)


            # Plot the data if telem (and tracestats) were defined
            if telem:
        
                tracestats['stt_' + tracedata.name] = tracedata

                # Define color to use
                tcolor = self.decplots['plots'][plotnum]['traces'][tracenum]\
                                      ['TCOLOR']
                tcolor = self._getcolor(tcolor=tcolor, tracenum=tracenum)        
        
                ax.plot(telem.plottimes, telem.plotdata, color=tcolor, 
                        label=tracedata.name, linewidth=1)  
                plt.draw()

    

        # Add Y label
        ax.set_ylabel(self.decplots['plots'][0]['PYLABEL'],
                      color=self.plotinfo['fgcolor'])

        # Configure y axis
        self._configureyaxis(ax, tracestats=tracestats)

        # Configure x axis
        self._configurexaxis(ax, plotnum)

        # Configure legend if data exists
        if tracestats:

            # Configure legend
            self._generatelegend(ax, plotnum)


        # Add a grid
        self._generategrid(ax)

        return tracestats




    def _plotlttaxis(self, plotnum):

        # Define attributes to keywords to simplify
        bgcolor = self.plotinfo['bgcolor']
        fgcolor = self.plotinfo['fgcolor']
        ax = self.axlist[plotnum]
        lttax = self.lttaxlist[plotnum]

        lttax.set_title('Long Term Trends', fontsize=12, color=fgcolor,
                        position=[0.5, 1.01])

        tracestats = {}
        for tracenum in self.decplots['plots'][plotnum]['traces'].keys():

            tstat = self.decplots['plots'][plotnum]['traces'][tracenum]['TSTAT']
            
            oneyearago = DateTime(self.time2).secs - self.ltttimespan
            oneyearago = DateTime(oneyearago).date
            telem, tracedata = self._getplotdata(plotnum, tracenum, 
                                                  time1=oneyearago,
                                                  fetchstat='daily')

            # Plot the data
            if telem:

                tracestats['ltt_' + tracedata.name] = tracedata
        
                # Define color to use
                tcolor = self.decplots['plots'][plotnum]['traces'][tracenum]\
                                      ['TCOLOR']
                tcolor = self._getcolor(tcolor=tcolor, tracenum=tracenum)      

                # When plotting the LTTs, disregard the plottimes and 
                # plotdata, instead plot the mins-to-maxes for each day,
                # unles this is a state based MSID or if a particular stat
                # was requested, in which case you want to use the plottimes 
                # and plotdata.
                if (telem.type == 'state') or tstat:
                    # Generate the trace data from the means
                    plottimes = telem.plottimes
                    plotdata = telem.plotdata
                else:
                    # Generate the trace data from the mins and maxes
                    plottimes = np.array(zip(telem.times, telem.times + 1))
                    plottimes = plottimes.flatten()     
                    plotdata = np.array(zip(telem.mins, telem.maxes))
                    plotdata = plotdata.flatten()

                lttax.plot(plottimes, plotdata, color=tcolor, 
                           label=tracedata.name, linewidth=1)  
            
            time1secs = DateTime(self.time1).secs
            time2secs = DateTime(self.time2).secs
            lttax.axvspan(time1secs, time2secs, color=[.6, .6, .6], alpha=0.2)

            t1 = DateTime(self.time2).secs - self.ltttimespan
            self._configurexaxis(lttax, plotnum, t1=t1, numticks=2)        
        
        #
        # This was moved to the figure level function, since it requires info 
        # from both the short term and long term plots/axes.
        #
        #if tracestats:
        #    
            # Configure y axes for both primary and ltt plots to be equal
            # self._configureyaxis(lttax, tracestats=tracestats)
            # self._configureyaxis(ax, tracestats=tracestats)

        self._generategrid(lttax)

        return tracestats





    def _plotbinaryaxis(self, plotnum):

        bgcolor = self.plotinfo['bgcolor']
        fgcolor = self.plotinfo['fgcolor']
        bax = self.baxlist[plotnum]
        ax = self.axlist[plotnum]

        msidkey = {}
        for tbtracenum in self.decplots['plots'][plotnum]['tbtraces'].keys():

            # Fetch telemetry and return name of MSID/Data
            tbtelem, tbtracedata = self._getplotdata(plotnum, tbtracenum, 
                                                     binaryplot=True)

            offval = tbtracenum * 2
            onval = offval + 1
            
            # First set all the values to stepped up off value
            bin = np.array([offval] * len(tbtelem.vals))

            # Then set all the on values to the stepped up on value
            bin[tbtelem.vals == 1] = onval  

            # Plot the data
            if tbtelem:
        
                # Define color to use
                tcolor = self.decplots['plots'][plotnum]['tbtraces']\
                                      [tbtracenum]['TCOLOR']
                tcolor = self._getcolor(tcolor=tcolor, tracenum=tbtracenum)        


                bax.plot(tbtelem.times, bin, color=tcolor, 
                         label=tbtracedata.name, linewidth=1)

                msidkey[tbtracenum] = tbtracedata.msid



        # Configure the y axis labels for the binary data (heater
        # states)

        bax.set_xticklabels('')
        bax.set_xlim(tbtracedata.mintime - 360, tbtracedata.maxtime + 360)

        # Configure binary trace y axis
        ymaxlim =  np.max(msidkey.keys()) * 2
        ytick = range(0, ymaxlim + 2, 2)
        ylab = [''] * len(ytick)
        for n in msidkey.keys():
            ylab[n] = msidkey[n]

        bax.set_ylim(-1, ymaxlim+3)
        bax.set_yticks(ytick)
        bax.set_yticklabels(ylab, fontsize=6, color=self.plotinfo['fgcolor'])
        bax.tick_params(axis='y', which='major', length=6, direction='out',
                        color=[0.8, 0.8, 0.8])

        # Reposition title since the reshaping of the figure affects
        # the y scaling of where the title is located.
        ax.set_title(self.decplots['plots'][plotnum]['PTITLE'], fontsize=12,
                     color=fgcolor, position=[0.5, 1.2])                




    def _plotfigure(self):
        
        # Define figure attributes to keywords to simplify
        bgcolor = self.plotinfo['bgcolor']
        fgcolor = self.plotinfo['fgcolor']
        figsize = (self.plotinfo['width'], self.plotinfo['height'])
        

        # Create figure
        plt.rc('axes', edgecolor=fgcolor)
        fig = plt.figure(figsize=figsize, facecolor=bgcolor)
        self.fig = fig


        #---------------------------------------------------------------------
        # Figure Annotations

        # Create title and subtitle
        pagetitle = fig.text(0.5, 0.98, self.decplots['DTITLE'], ha="center",
                             va="center", size=14, color=fgcolor)
        pagesubtitle = fig.text(0.5, 0.96, self.decplots['DSUBTITLE'], 
                                ha="center", va="center", size=12, 
                                color=fgcolor)

        # Create source file info (for upper lefthand corner)
        decname = 'Filename: %s\n'%os.path.basename(self.decfile)
        timeperiod = 'Date Range: %s to \n            %s\n'%(self.time1,
                                                            self.time2)
        generated = 'Generated: %s\n'%DateTime().date
        datasource = 'From: Engineering Telemetry Archive'

        sourceinfo = fig.text(0.01, 0.99, decname + timeperiod + generated +
                              datasource, ha="left", va="top",
                              size=self.plotinfo['stampfontsize'],
                              family='monospace', color=fgcolor)


        #---------------------------------------------------------------------
        # Axis Handle Placeholders

        # Create empty axes lists for the primary axis, the binary sub
        # axes (in case they are required), and LTT axes.
        self.axlist = [None for n in range(self.decplots['numplots'])]

        self.baxlist = [None for n in range(self.decplots['numplots'])]

        if self.plotltt:
            self.lttaxlist = [None for n in range(self.decplots['numplots'])]
        

        #---------------------------------------------------------------------
        # Main loop for generating all plots 

        for plotnum in np.sort(self.decplots['plots'].keys()):

            print('plot %d = %s'%(plotnum,
                                  self.decplots['plots'][plotnum]['PTITLE']))
            

            # Create Main Axes
            plotloc = self.plotinfo['location'][plotnum]
            self.axlist[plotnum] = fig.add_axes(plotloc, axisbg=bgcolor)


            # Create Primary Plots
            tracestats = self._plotaxis(plotnum)


            # Configure y axes for primary plot
            self._configureyaxis(self.axlist[plotnum], tracestats=tracestats)


            if self.plotltt:
                lttplotloc = self.plotinfo['lttslocation'][plotnum]
                self.lttaxlist[plotnum] = fig.add_axes(lttplotloc, axisbg=bgcolor)
                
                ltttracestats = self._plotlttaxis(plotnum)

                if ltttracestats:
                    # Use the LTT stats, this will be used to scale the y axes.
                    # It is assumed that the LTT plot will include the primary
                    # plot data, so the LTT stats (extrema) can be used 
                    # instead.
                    tracestats.update(ltttracestats)

                # Configure y axes for both primary and ltt plots to be equal
                self._configureyaxis(self.axlist[plotnum], 
                                     tracestats=tracestats)
                self._configureyaxis(self.lttaxlist[plotnum], 
                                     tracestats=tracestats)


            # Create binary axes if defined (i.e. if it isn't None)
            if self.plotinfo['binarylocation'][plotnum]:
                binaryplotloc = self.plotinfo['binarylocation'][plotnum]
                bincolor = axisbg=self.plotinfo['bgcolor']
                self.baxlist[plotnum] = fig.add_axes(binaryplotloc, 
                                                     axisbg=bincolor)

                self._plotbinaryaxis(plotnum)



            # Write the stats text for the current plot to the right of the 
            # plot.
            self._writestats(tracestats, plotnum)

            
        # Plot date labels in figure coordinate system
        #
        # plotloc should correspond to the bottom plot
        if self.plotinfo['binarylocation'][plotnum]:
            self._plotdatelabels(binaryplotloc)
        else:
            self._plotdatelabels(plotloc)

        if self.plotltt:
            lttstart = DateTime(self.time2).secs - self.ltttimespan
            self._plotdatelabels(plotloc, t1=lttstart, t2=None, numticks=2, 
                                 lttplot=True)

        # # Format the date axis
        # fig.autofmt_xdate(rotation=0, ha='center')


        #---------------------------------------------------------------------
        # Save the plot to a file

        basename = os.path.basename(self.decfile).split('.')[0]

        t1 = DateTime(self.time1).greta
        t2 = DateTime(self.time2).greta
        
        filename = os.getcwd()+'/' + basename + '_' +t1 + '-' + t2 + '.png'

        fig.savefig(filename, facecolor=self.plotinfo['bgcolor'],
                    edgecolor=self.plotinfo['bgcolor'])

        print('Saved plot to %s'%filename)

               
        # Add ability to manually adjust y axis

