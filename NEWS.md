# eegr 0.3-2
- date: 20220511
- added trd objects and methods
- added mutitaper spectral estimation: mtspec(), mtchan(), mtfft()
- added cleanline(), removetrend(), ssinterp(), unusablesensors(),
        noisysensors(), cprep(), eeginterp(), reref() functions
- added 'ref'' attribute to ctd object
- used inherits() instead of class()

---

# eegr 0.3-1
- date: 20211212
- bugfix in readEDF()

---

# eegr 0.3-0
- date: 20210602
- use Akima spherical splines in topoplot()
- added changefs(), selectdata(), filterdata(), detrenddata(), 
- added sample dataset EEGdata
- topoplot: check for missing values in x (produced error in predict.smooth.sspline)
- topoplot: repaired bug that plotted no sensorlocs
- ctd: functions fs(), ns(), npts() return numeric value, not num [1(1d)]
- consistent use of GPLv3
- replaced signal by gsignal

---

# eegr 0.1.1
- date: 20191210
- bugfix in readEDF(), error reading empty annotation channel

---

# eegr 0.1.0
- date: 20190918
- Added functions:
- ctd - object and methods
- sensorlocs - object and methods
- readSensorLocs
- eeglocations data set 
- readEDF
- topoplot

---

# eegr 0.0.5
- date: 
- create eegr_utils.R
- move function unfactor (sensorlocs.default) to eegr_utils.R
- reset plotting parameters in plot.sensorlocs
- define functions taper, 
- add tests for taper, spectrum

---

# eegr 0.0.4
- date: 20160706
- bugfix get.sensorlocs.labels
- define method plot.ctd.data

---

# eegr 0.0.2
- date: 20160626
- add EEGlocations data
- define object and methods for sensorlocs
  = default, plot (requires plotrix library), get.sensorlocs.labels
- define function readsensors

---

# eegr 0.0.1
- date: 20160609
- define object and methods for ctd.data
  = default, print, summary, print.summary
- define function readbdf (requires mmap library)

