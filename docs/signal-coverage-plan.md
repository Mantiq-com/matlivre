# Signal Processing Toolbox: coverage plan

How MatLivre grows to cover the MATLAB Signal Processing Toolbox. The function inventory is
`tools/coverage/signal.tsv` (category, name), taken from the toolbox's function reference
(mathworks.com/help/signal/referencelist.html, R2026b listing). Functions the reference marks
"Since R2024a" or later are left out, because MatLivre targets R2023b.

## Where things stand

| | Count |
|---|---:|
| Functions in scope (R2023b-era, from `signal.tsv`) | 333 |
| Already in MatLivre before this plan (`xcorr`, `xcov`, `resample`, `rms`, `sinc`, `detrend`, `movmedian`, `movmad`, `isoutlier`, `filloutliers`, `smoothdata`, `islocalmax`, `paddata`, `trimdata`, `resize`, `freqspace`, `ss2tf`) | 17 |
| Implemented by this plan (`crates/livre-lang/src/signal/`) | 237 |
| Still to implement and verifiable here | 65 |
| Not verifiable: requires Wavelet/Audio Toolbox, not installed (see below) | 14 |

The 237 are checked against MATLAB: see "How it is tested".

## Not verifiable: requires Wavelet/Audio Toolbox, not installed

The MATLAB install used for recording has neither the Wavelet nor the Audio Toolbox, so these cannot be
recorded and are not implemented: `cwt`, `modwt`, `modwpt`, `imodwt`, `wsst`, `wcoherence`, `tqwt`, `cqt`,
`ewt`, `spectralCrest`, `spectralEntropy`, `spectralFlatness`, `spectralKurtosis`, `spectralSkewness`.
They stay in the gap count; the `dwt`/`wavedec` family is likewise unrecordable.

## How it is tested

`crates/livre-lang/tests/all/signal_matlab.rs` runs `data/signal_cases.m` (902 calls covering valid
use, option variants, orientation, edge cases and error messages) and compares it with
`data/signal_matlab.txt`, recorded from MATLAB R2023b Update 3 with Signal Processing Toolbox 23.2.
Classes, sizes and error messages must match exactly and numbers to a relative 1e-9 (1e-6 for single
precision). 881 of the 902 cases match; the other 21 are listed in
`data/signal_known_deviations.tsv` with the reason, and the test fails if one of them starts to
match, so the list only shrinks.

To add a function: add its calls to `tools/signal/cases.txt` (each line assigns `r`, a number array
or a cell array of them), run `python tools/signal/gen_cases.py` to rebuild `signal_cases.m`, run
that in MATLAB (`matlab -batch "run('signal_cases.m')"`) and keep its result lines as
`signal_matlab.txt` (the generator's docstring has the commands). The same script runs here; any
difference is either a bug to fix or an entry for the known-deviations list.

Running MATLAB over the 76 functions found most of the conventions the reference leaves open
(window length 0 and 1, row or column orientation, zero padding of transformed filters, the order
of `chirp`'s quadratic shapes, automatic lengths for `impz` and `freqs`, the exact wording of
validation messages) and several outright mistakes in the first, documentation-only
implementation. Those are now fixed and covered by the recorded cases.

## Out of scope

Following the existing scope rules for a live-script engine:

- Apps and interactive tools (Signal Analyzer, Signal Labeler, Signal Feature Extractor, Filter
  Designer, Filter Analyzer, Window Designer, WVTool, Wavelet and Multiresolution analyzers,
  EDF File Analyzer) and the Live Editor tasks (Design Filter, Create Plot).
- Deep-learning layers and datastores: `signalDatastore`, `cwtLayer`, `stftLayer`,
  `modwtLayer`, `dlcwt`, `dlstft`, `dlmodwt`, `deepSignalAnomalyDetector`.
- `filt2block` (Simulink) and code generation.
- The Slepian database functions `dpssclear`, `dpssdir`, `dpssload`, `dpsssave`
  (`dpss` itself is in scope).
- Signal labeling objects and region-of-interest functions (15 functions, object-heavy).
- EDF file I/O (`edfheader`, `edfinfo`, `edfread`, `edfwrite`): a file-format task.

That removes about 32 of the 333 names from the practical target; they stay in `signal.tsv`
so the coverage page can list them as out of scope.

## Phases

Each phase is a set of modules under `crates/livre-lang/src/signal/`, wired like the other
builtin modules (a `NAMES` list per module, `signal::call` as the router). New functions get
cases in `signal_cases.m` and MATLAB's output recorded as above.

1. **Windows, decibel conversions, basic waveforms. Done**: `rectwin`, `bartlett`, `triang`,
   `hamming`, `hann`, `blackman`, `blackmanharris`, `nuttallwin`, `flattopwin`, `kaiser`,
   `tukeywin`, `gausswin`, `parzenwin`, `bohmanwin`, `barthannwin`, `chebwin`; `mag2db`, `db2mag`,
   `pow2db`, `db2pow`; `square`, `sawtooth`, `rectpuls`, `diric`, `chirp`, `gmonopuls`.
2. **More waveforms and windows**. Done: `tripuls`, `gauspuls`, `buffer`, `shiftdata`, `unshiftdata`, `uencode`, `udecode`, `taylorwin`,
   `pulstran`, `vco`, `modulate` (`pwm` and `ppm` are not implemented), `dpss` (direct method only), `marcumq`, `demod` (`pwm` and `ppm` are not implemented).
3. **Polynomial and filter conversions**. Done: `lp2lp`, `lp2hp`, `lp2bp`, `lp2bs` (transfer
   function form), `bilinear` (transfer function and zero-pole-gain forms, with pre-warping),
   `polystab`, `polyscale`, `eqtflength`, `sos2tf`, `tf2zp`, `tf2zpk`, `zp2tf`, `sos2zp`, `zp2sos`, `tf2sos`, `tf2ss`, `zp2ss`, `ss2zp`, `sos2ss`, `ss2sos`, `residuez`, `latc2tf`, `tf2latc`, `impinvar`, `cell2sos`, `sos2cell`.
4. **Analog prototypes and IIR design**. Done: `buttap`, `cheb1ap`, `cheb2ap`, `butter`,
   `cheby1`, `cheby2`, `ellipap`, `besselap`, `ellip`, `besself`, `buttord`, `cheb1ord`, `cheb2ord`, `ellipord` (digital and analog, all four band
   types, transfer function, zero-pole-gain and state-space outputs), `maxflat` (without the second-order-section
   outputs). To do: `yulewalk`.
5. **FIR design and filter analysis**. Done: `fir1`, `freqz`, `freqs`, `grpdelay`, `phasez`, `phasedelay`, `zerophase`,
   `isstable`, `impz`, `stepz`, `sgolay`, `fir2`, `firls`, `kaiserord`, `firpm`, `firpmord`, `gaussdesign`, `rcosdesign`, `intfilt`, `impzlength`, `filternorm`, `filtord`,
   `firtype`, `isminphase`, `ismaxphase`, `islinphase`, `isallpass`, `fircls`, `fircls1` (without the WT and
   WP/WS/K options), `zplane` (zeros, poles, multiplicity labels and axes of the plot; not the `digitalFilter` input).
   `cfirpm` (amplitude vectors, your own response functions, and the predefined `@lowpass`, `@highpass`, `@bandpass`, `@bandstop`, `@multiband` without a delay, `@differentiator`, `@hilbfilt` and `@invsinc`).
   To do: `designfilt`, `digitalFilter`.
6. **Filtering and multirate**. Done: `decimate`, `interp`, `latcfilt`, `filtfilt` (transfer function form), `sosfilt`, `filtic`,
   `medfilt1`, `downsample`, `upsample`, `fftfilt`, `upfirdn`, `sgolayfilt`, `hampel`. To do:
   `lowpass`, `highpass`, `bandpass`, `bandstop`, `fillgaps`;
   `filtfilt` with second-order sections.
7. **Spectral estimation and transforms**. Done: `periodogram`, `pwelch`, `mscohere`, `cpsd`, `tfestimate`, `dct`,
   `idct`, `dftmtx`, `goertzel`, `hilbert`, `cconv`, `convmtx`, `bitrevorder`, `czt`, `fwht`, `ifwht`, `cceps`, `rceps`, `icceps`,
   `digitrevorder`, `corrmtx`, `pyulear`, `pburg`, `pcov`, `pmcov`, `xcorr2`, `pmusic`, `peig`, `rootmusic`, `rooteig`, `pmtm`, `plomb`, `poctave` (power and PSD smoothing; not `Weighting` or the spectrogram type).
   `pspectrum` (power and spectrogram types), `modalfit`, `modalsd` (see their notes below).
8. **Measurements and modeling**. Done: `peak2peak`, `peak2rms`, `rssq`, `findpeaks`, `enbw`, `bandpower`, `meanfreq`, `medfreq`,
   `finddelay`, `seqperiod`, `cusum`, `dtw`, `edr`, `alignsignals`, `obw`, `powerbw`, `statelevels`, `midcross`, `risetime`, `falltime`,
   `pulsewidth`, `pulseperiod`, `pulsesep`, `dutycycle`, `overshoot`, `undershoot`, `settlingtime`, `slewrate`, `zerocrossrate`, `findchangepts`, `envelope`, `invfreqz`, `invfreqs` (non-iterative fit), `stmcb`, `findsignal` (fixed time alignment only). To do:
   `snr`, `sinad`, `thd`, `sfdr`, `toi`;
   Done: `levinson`, `ac2rc`, `rc2ac`, `rc2poly`, `poly2rc`, `ac2poly`, `poly2ac`, `rc2is`, `is2rc`,
   `rc2lar`, `lar2rc`, `poly2lsf`, `lsf2poly`, `aryule`, `arburg`, `arcov`, `armcov`,
   `rlevinson`, `schurrc`, `prony`, `lpc`.
9. **Time-frequency**. Done: `spectrogram`, `xspectrogram`, `ifsst`, `tfridge`, `iscola`, `stft`, `istft`, `wvd` (not `smoothedPseudo`), `instfreq` (both methods, `tfd,fd,td` input), `instbw`. Not recovered after one round: `xwvd`. `vmd` (explicit ascending `CentralFrequencies` exact; default initialisation not recovered). To do: `xspectrogram`,
   `hht`, `emd`, `vmd`, `ewt`, `fsst`, `kurtogram`,
   spectral descriptors (`spectralCrest`, `spectralEntropy`, `spectralFlatness`,
   `spectralKurtosis`, `spectralSkewness`) and the wavelet functions (`cwt`, `modwt`,
   `modwpt`, `wsst`, `wcoherence`, `tqwt`, `cqt`): the last two groups are not verifiable here.
10. **Vibration analysis**. Done: `rainflow`, `envspectrum`, `tsa`, `tachorpm` (FitType 'linear' only). To do: `rpmtrack`, `orderspectrum`, `ordertrack`,
    `orderwaveform`, `rpmordermap`, `rpmfreqmap`, `modalfrf`,
    `modalfit`, `modalsd`.

## Known deviations from MATLAB

These are the calls in `signal_known_deviations.tsv`, plus the limits that the recorded cases do
not exercise:

- **Order of zeros and poles** of `cheby2` designs. `butter`, `cheby1` and the poles of `cheby2`
  follow MATLAB's order, because the design goes through the same state-space path (cascade
  realisation `zp2ss`, frequency transformation, bilinear transform) and takes the eigenvalue
  order of the state matrix from a port of LAPACK's `dgebal`/`dgehd2`/`dlahqr`
  (`signal/lapack_eig.rs`). The zeros of `cheby2` are the eigenvalues of A - B C / D; in about one
  case in twenty MATLAB's order follows rounding errors in the entries of that matrix, which
  differ here (five recorded cases).
- **Degenerate designs**: `cheby1` with zero ripple and `cheby2` with zero attenuation.
- **Single precision**: results are computed in double precision and rounded to single; MATLAB's
  own single arithmetic can differ in the last digits (the recorded cases allow 1e-6).
- **N-D arrays** and operating along dimension 3 (`dct`, `hilbert`).
- **`spectrogram`**: the centroid outputs `fc`, `tc`, and a frequency vector in place of `nfft`.
- **`pwelch`** with a non-integer overlap; **`freqz`** with a zero denominator;
  **`filtic`** with a zero leading denominator coefficient; **`gauspuls`** with complex time.
- **`grpdelay`** near zeros on the unit circle: MATLAB's own value loses accuracy there (3e-7
  off at the last point of a 5th-order Butterworth); MatLivre's is within 1e-8 of the exact one.
- **`dct`/`idct`** sum the definition directly, O(N^2) per column; fine for a few thousand samples,
  slow beyond that.
- **`fftfilt`** and **`upfirdn`** filter directly instead of with overlap-add FFTs or a polyphase
  structure; same results up to rounding.
- **`filtfilt`** takes a transfer function only (no second-order sections).
- **`fircls`/`fircls1`**: MatLivre returns the exact constrained optimum; MATLAB stops when the bounds
  are violated by under about 5e-9, so its coefficients differ by 1e-7 to 2e-6 (every numeric
  recorded call is a listed deviation for this reason). The WT and WP/WS/K options of `fircls1` and
  zero lower bounds in `fircls` are not implemented.
  One extra round tried a stopping tolerance of 5e-9 on the exact active-set iteration (the exchange
  iteration of the 1996 paper is not reproducible step for step from the published description);
  no recorded call moved closer, so the exact optimum stays.
- **`pmusic`/`peig`**: the noise eigenvectors (third output) agree with MATLAB up to the sign of each
  column (or the basis of a nearly degenerate noise subspace); the pseudospectrum, frequencies and
  eigenvalues agree. MATLAB's data matrix holds each segment in reverse order, and for an explicit
  window length the last segment is zero padded and the eigenvalues are normalised by one less than
  the segment count; `peig` weights use a different normalisation again (see `subspace.rs`).
- **Deferred queue** (after the research items): `cfirpm` predefined responses (`@lowpass`, `@multiband`, ...); the `line` gaps moved to `docs/known-limits.md`.
- **Test time**: the cfirpm designs of order 16 and above (about 9 s in a debug build) run only in the
  ignored test `signal_heavy_cases_match_matlab` (`cargo test -p livre-lang --test all -- --ignored`);
  the regular `signal_functions_match_matlab` skips them.
- **Research queue**: `tachorpm` 'smooth' fit (least-squares spline with FitPoints breakpoints); `vmd` initialisation ('peaks' start picks the largest spectral peaks of the signal, order and exact picking unknown), update order for non-ascending starts, odd-length mirror; `emd` and `vmd` (one `emd` round: spline envelopes through the extrema with or without the end points, Cauchy sifting stop at 0.2, did not reproduce the first IMF; `hht` depends on `emd`); `xwvd` and `wvd(...,'smoothedPseudo')` (plain `wvd`: products x(a) y*(b) with a+b = 2t and lags |m| <= N/2 at N frequencies k/(2N) cycles per sample, real input made analytic from the record zero padded to 2N; the cross distribution does not follow that rule); `pspectrum` (power type only first): the output is always 4096 points from 0 to Nyquist, equal to twice the Kaiser (beta 40(1-leakage)) windowed modified periodogram of the whole record for short signals, but it deviates from the exact DTFT of that window by up to 1e-4 between some frequencies (an interpolation from a coarser grid?) and the Welch regime for long records (more than about 20000 samples) was not probed; `poctave` low-band filters (one round: decimate-then-filter and an independent bilinear design did not match); `rootmusic`/`rooteig` powers (systematic ~1e-7 relative; one round tried: mean noise eigenvalue as the noise variance, p/2 or p signal eigenvalues, exact or least squares; the zero-frequency root case gave P=[-0.24, 44.6] with a 2x2 solve, not MATLAB's [0, 34.29]).
- **`rootmusic`/`rooteig`**: frequencies agree to 1e-12, powers to about 1e-7 relative (MATLAB's
  power system was not recovered); a real correlation matrix with a root at zero frequency returns
  one power where MATLAB returns two.
- **`tachorpm`**: pulses are the means of consecutive edge pairs (edge times = middles of the risetime and
  falltime crossings between the state levels), speeds 60 / (interval * PulsesPerRev) at interval
  midpoints, linearly interpolated (and extrapolated) onto (0:N-1)/OutputFs. The default FitType 'smooth'
  is not implemented (research queue).
- **`vmd`**: ADMM on the one-sided spectrum of the signal mirrored by half its length on each side
  (mode update 1 / (1 + alpha (f - f_k)^2) with f in cycles per sample, Lagrange update with the
  mean-square change over the bins, stop when the absolute improvement (sum / bins) < AbsoluteTolerance
  and the relative one < RelativeTolerance, modes returned by decreasing central frequency, two-sided
  reconstruction with DC and Nyquist counted once). Recovered exactly for explicit ascending
  `CentralFrequencies`; the default 'peaks', 'grid' and 'random' starts, descending starts, odd lengths
  and single precision are listed deviations (research queue).
- **`tsa`**: cycles run between pulses (a scalar `tp` is an interval; pulses reach the record duration
  N/fs); one rotation is the shortest spell of PulsesPerRotation pulses, `M = NumRotations * floor(rotation * fs + 0.5 - 1e-9)`,
  t = (0:M-1)/M * period, windows of NumRotations rotations do not overlap, each rotation is interpolated onto
  M/NumRotations points (the 'fft' method cuts whole-sample segments and averages their spectra), and
  ResampleFactor first upsamples the signal by Fourier interpolation. Single precision differs at 1e-6.
- **`poctave`**: the ANSI S1.11 center frequencies; bands that overlap the limits (the band holding
  Nyquist is kept while its lower edge is under 0.9 of the upper limit and is a highpass of order
  FilterOrder); each band is a Butterworth bandpass of order FilterOrder (the SOS form) and the
  power is the mean square of its output. Bands with a center frequency below about 1/250 of the
  sample rate do not match (see the research queue), and `Weighting` and the spectrogram type are
  not implemented. PSD smoothing (rectangle rule, edge bins at 0 and Nyquist counted twice) matches.
- **`plomb`**: MATLAB uses the Press-Rybicki fast algorithm (4-point extirpolation, a power-of-two
  FFT, period ofac*N*ts), so MatLivre does too; a given frequency vector is evaluated exactly. The
  spectrum is scaled by 2 var ts (psd; ts from the whole record), 2 var / N (power) or 1; the
  detection thresholds use N*fmax/fNyquist independent frequencies. Differences: single precision and
  two-to-five-sample records (rounding).
- **`pmtm`**: the adaptive combination (Thomson weights from an initial two-taper average) stops
  when the summed change of the spectrum over the FFT grid falls below 5e-4 times the signal power
  (the rule that fits every recorded call); the confidence bounds use 2k degrees of freedom for all
  three combinations; with `ConfidenceLevel` the outputs come as `[pxx,w,pxxc]`. Not matched: nothing
  recorded; sine tapers, matrix input, complex input and custom frequencies are covered.
- **`cfirpm`**: MATLAB's result is the exact Chebyshev optimum on its frequency grid (the band edges
  plus the multiples of 2/L at least half a step from both edges; L = 2^nextpow2(0.8 lgrid (n+1))
  for the default lgrid of 25, which other lgrid values only approximate; for SYM='none' bands from -1
  to 1 the last two points become 1 - 1/(lgrid (n+1)); `real` drops the point one step below 1).
  MatLivre finds that optimum with a barrier method. Listed deviations: positive-frequency-only
  complex designs (the optimum is not unique), results MATLAB stops short of the optimum by about
  1e-6, negative band weights, complex amplitude values, and the `res` extremal lists. The
  predefined responses (`@lowpass`, `@bandpass`, ...) raise a "not implemented" error. Odd-symmetric
  designs from amplitude vectors return a zero filter, as MATLAB does.
- **`zplane`**: the unit circle has 100 points, more as a root nears the circle (round(1/distance),
  at most 50000); the axis limits are the 10 percent padded extents of everything drawn, widened to
  the axes box ratio. Not matched: single precision data (MatLivre's line objects store doubles), the
  limits left by a second call with hold on, and axes made by subplot (their Position is not tracked).
  `zplane(d)` for a `digitalFilter` is not implemented.
- **`maxflat`**: the second-order-section outputs (`sos`, `g`) are not implemented (MATLAB's section
  order and pairing were not reproduced); orders around 30 and the degenerate Wn = 0.5 designs
  differ at 1e-6 to 1e-8 because the flatness system is ill-conditioned.
- **`impz`/`stepz`** choose their automatic length from the poles with the rule inferred from
  MATLAB's output (floor of ln(5e-5)/ln(radius) for stable poles, ln(1e6)/ln(radius) for unstable
  ones, five periods on the unit circle); unusual mixtures of pole radii may differ by a sample.
- **`freqs`** chooses its automatic range from the root magnitudes (a root at the origin counts
  as 1) by a rule inferred from MATLAB's output. MATLAB also refines the grid around resonances
  (second-order or higher denominators), which MatLivre does not: it spaces the points evenly.
- **`hilbert`** ignores the imaginary part of complex input, with MATLAB's warning.
- **Messages**: the wording follows MATLAB's for the cases recorded; other validation failures
  fall back to generic text.

## Per-phase checklist

From `AGENTS.md`:

1. Implement in `src/signal/<area>.rs`; add the area's `NAMES` to the combined list in
   `builtins.rs` and to the router in `signal/mod.rs`.
2. Add calls to `tests/all/data/signal_cases.m`, record MATLAB's output into
   `signal_matlab.txt`, and fix or list every difference.
3. Build and run tests one at a time (shared machine; no other `cargo` running).
4. Commit the implementation, run `bash tools/coverage/update.sh`, commit "Refresh the
   coverage checklist".
5. Update `README.md`: function count and the feature lists.

## Tooling still to do

- Feed `tools/coverage/signal.tsv` into `update.sh` / `coverage.pl` so
  `docs/matlab-coverage.md` gains a "Signal Processing Toolbox" section. `documented.tsv` is
  generated from the installation's `helpfuncbycat.xml`; the toolbox's reference list is now
  available in the MATLAB installation, so it can be generated the same way.
- A test that every name in each module's `NAMES` is dispatched is covered by the recorded cases
  for the implemented functions; new functions need their cases.

### modalfrf (done, partial)
Single fixed input with one or more outputs, `Measurement` rovinginput/rovingoutput, H1/H2 exact, Sensor acc/vel/dis, coherence output. Not recovered: Hv estimator (MATLAB's formula differs from the total-least-squares form and from sqrt(H1*H2)), multiple fixed inputs (MIMO), subspace, non-integer window/overlap. modalfit/modalsd still to do.

### rpmfreqmap, rpmordermap (done)
Both exact against MATLAB (single-precision `rpmordermap` listed as a deviation). `rpmfreqmap`: symmetric window of the shortest length whose resolution `fs*sum(w^2)/sum(w)^2` does not exceed the target (default `fs/128`), `nfft = max(256, L)`, hop `floor(L*(1-ov))`, rms/peak/power scaling with DC and Nyquist bins special, rpm interpolated at the segment centre. `rpmordermap`: phase by `cumtrapz` of the 15x resampled rpm, linear interpolation onto `fp = 8*Omax` samples per revolution, `res = fp/256`, window length solved the same way, `noverlap = min(ceil(ov*L), L-1)`, `floor(L/8)+2` order rows. This also added the numeric `resample(x,p,q,n,beta)` (Kaiser-windowed sinc, `p*fir1`-style DC gain) in front of the timeseries `resample`.

### Not recovered (one round each)
- `orderwaveform`, `ordertrack`, `rpmtrack`: Vold-Kalman filter. A plain first/second-order structure-equation least-squares solution with one weight `r` leaves a 10-15 % residual against MATLAB for every `r`, so the weighting/bandwidth mapping was not recovered; `rpmtrack` also needs the coarse STFT ridge step.
- `orderspectrum`: one round earlier (order axis of 122 points, `dO = 0.24819*fs/rpm`); still not recovered.

### pspectrum (done: power and spectrogram types)
Exact against MATLAB except the cases in the deviations file. Findings: Kaiser window with beta = 40*(1-Leakage); the equivalent noise bandwidth MATLAB uses is not the discrete one (table read from the limits in its error message, 0.01 steps, cubic interpolation between); output grid = linspace over the limits with 4096 points (1024 for spectrograms) or `round(8*fnyq/FR)+1` points when `FrequencyResolution` is given; spectra are DTFTs on that grid (folded FFT), one-sided values doubled at every bin including DC. Power: a single modified periodogram when the target resolution equals `fs*ENBW/N`, else Welch with `L = round(fs*ENBW/(8*df))`, hop `L/(2*ENBW-1)`, `K = ceil((N-L)/hop)+1` segments at `round(linspace(0,N-L,K))`; the default floor on the resolution is `fs/1024`. Spectrogram: `L = ceil(N/d)` with d = 2 (N<64), 8, 16 (N>=2048), 32 (4096), 64 (8192), 128 (16384); hop `ceil(L*(1-ov))`; the last segment is shortened and gets its own window of that length. Not implemented: `persistence`, `Reassign`, nonuniform time vectors, timetable input.

### modalfit, modalsd (done: LSCE, single FRF)
`modalfit` (fn, dr, ofrf exact; `ms` approximate): impulse response = real inverse FFT of the one-sided FRF (length `2(p-1)`, DC and Nyquist real), covariance-method Prony of order `2*mnum` on the first `22*mnum` samples (`mldivide`), poles `log(roots)*fs`; natural frequencies from the poles with positive imaginary part of `z` (real roots are dropped and the vectors padded with NaN to `mnum`); residues from a frequency-domain least-squares fit of `1/(jw-s)` over all `2*mnum` poles; `ofrf` is the sum over the identified modes. `ms` = square root of those residues agrees with MATLAB to 1e-5 when the fit is good and not when it is poor (MATLAB's mode shape normalisation was not recovered; also `PhysFreq` re-fits them). Not implemented: `FreqRange`, `lsrf`, `pp`, several FRFs (p-by-m-by-n). `modalsd`: row `i` of the `MaxModes`-square matrix holds the order-`i` poles that are stable against the nearest pole of order `i+1` (relative frequency change under `SCriteria(1)`, absolute damping change under `SCriteria(2)`), NaN otherwise, zeros above the diagonal; `MaxModes` is capped at `floor(2(p-1)/22)` with a warning. MATLAB also opens a (hidden) figure when `modalsd` is called with an output; that is not reproduced.

### Research queue, one more dense round (not recovered)
- `emd`: spline envelopes through the extrema (plain, end points added, linear end extrapolation, mirrored extension) do not reproduce the first sift; the mean envelope agrees to 1e-4 in the interior, so the end handling is what is missing.
- `snr` family: for a clean sinusoid the signal power is exactly `sum(P(lo:hi))/ENBW` over the monotone mainlobe of a Kaiser(beta 38) periodogram (power-scaled), but the noise estimate (median level, `min(level, P)`) was not recovered.
- `fillgaps`: forward and backward Burg predictions weighted linearly across the gap reproduce MATLAB to 6e-3; the exact fit (which samples feed each model) was not recovered.
- `filtfilt` with second-order sections: interior samples equal the transfer-function version to 1e-12, the edge samples differ by up to 2.5e-2 from every padding/initial-condition variant tried.
- Vold-Kalman (`orderwaveform`, `ordertrack`, `rpmtrack`): the coupled real least-squares form (data weight 1, structure weight `r`) with first and second order differences was tried; best fits leave 15-20 % error for any `r`.
- The lowpass family needs `designfilt` (deferred).

### ifsst (done), fsst (not recovered)
`ifsst`: `x = (2*Re(sum of the used rows) - Re(DC row) - Re(Nyquist row, even window length)) / (L*g(L/2))` for a one-sided matrix (`nr = floor(L/2)+1`), `sum / (L*g(L/2))` for a two-sided one (`nr = L`), window default `kaiser(256,10)`; a frequency range keeps the rows from the first `f >= lo` through the first `f >= hi`; each ridge column keeps `NumFrequencyBins` (4) rows each side of the ridge index and the columns are independent (nothing is removed between them). `fsst` was taken to the point where 91 % of the cells agree: STFT with the symmetric `kaiser(L,10)` window, hop 1, zero padding, segment starting L/2 before the sample, phase reference at index L/2; frequency estimate `eta + Im(V_g'/V_g)/(2 pi)` with the negated analytic window derivative; all `L` two-sided bins squeezed to `round(L*estimate)` (half away from zero) and the one-sided rows kept. The cells whose estimate is near a rounding boundary (about 9 % of the cells, small-energy ones) still differ, so `fsst` is not shipped.

### tfridge (done), kurtogram (not recovered)
`tfridge`: penalised dynamic programme on `-ln|tfm|` (jump cost `penalty*(j-k)^2`, first index wins ties); with penalty 0 the per-column maximum; further ridges zero `NumFrequencyBins` bins either side of the previous path; `fridge = f(iridge)` keeps the class of `f`, `iridge` and `lridge` follow the class of `tfm`; no error for `NumFrequencyBins` above a quarter of the rows. `kurtogram`: output shape is `2*level` rows by `3*2^level` columns (default level 4 for 1024 samples), `w = [2 4 6 8 12 16 ...]`; the filter bank of the fast kurtogram (Antoni and Randall) was not reproduced, so it is not implemented.

### orderspectrum (done); research queue after one more round
`orderspectrum`: the order-RPM map of `rpmordermap` (rms, default window) averaged over time: `sqrt(mean(map.^2, 2))` for rms and peak maps, `mean(map, 2)` for a power map; the map syntax gives the same numbers as the signal syntax. Order axis = the map's.
Not recovered after one more dense round: `wvd(...,'smoothedPseudo')` and `xwvd` (default grids are the plain ones, 2*ceil(N/2) frequencies by 4*ceil(N/2) times; a time-window/frequency-window smoothing of the lag products with Kaiser(beta 20) windows of the documented default lengths did not reproduce the values), `vmd` 'peaks' (starting at the three largest `findpeaks` maxima of the one-sided |FFT| gives imfs within 3e-6 of MATLAB's but not exactly; the true starts are some other function of the spectrum), `rootmusic` powers (the frequencies match to 1e-12, the powers differ at 1e-7: the amplitude system is ill-conditioned, so a different solver is enough to move them).

### cfirpm predefined responses
`@lowpass`, `@highpass`, `@bandpass`, `@bandstop` and `@multiband` (no delay offset) are exactly the amplitude vectors `[1 1 0 0...]`, `[0 0 1 1]`, `[0 0 1 1 0 0]`, `[1 1 0 0 1 1]` and `a`; with negative band edges `@bandpass`/`@bandstop` swap their bands; `@highpass`/`@bandstop` raise an odd order by one. `@differentiator` is `j*pi*fs*f` with weight `w/|f|`, `@hilbfilt` is `-j*sign(f)`, both odd-symmetric (the coefficients are the antisymmetric half, `delta` is purely imaginary), `@invsinc` is `pi*a*f/sin(pi*a*f)` in the first band and 0 elsewhere (4 or 6 band edges). Not implemented (listed deviations): a group-delay offset `D` (the complex design is not unique and MATLAB's solution differs), `@allpass`, and `@hilbfilt` with negative edges.

### poctave low bands (one more round)
Only the band whose center is below fs/250 differs (band 3.98 Hz at fs 1000: 4.30e-5 against MATLAB's 4.02e-5; the next band agrees to 3e-10 and the others to 1e-12). Decimating by 2 to 16 before the filter does not reproduce MATLAB's value either. The sixth-order bandpass has its poles within 0.0035 of the unit circle, so coefficient rounding alone moves its output by percent; matching MATLAB would need its exact sequence of floating-point operations in the design, which cannot be recovered from outputs.

## Design note: designfilt and digitalFilter (not started)
**Object.** A `digitalFilter` is an `Object` value of a built-in class (as the containers and the ode objects are), with read-only properties in MATLAB's order: `Coefficients` (a numerator row vector for FIR; an L-by-6 SOS matrix for IIR designs), then the specifications `FrequencyResponse`, `ImpulseResponse` ('fir'/'iir'), `SampleRate` (default 2, normalised frequency), `FilterOrder`, the frequency properties of the design and `DesignMethod`. The display is the one recorded from MATLAB: a header `d = ` then ` digitalFilter with properties:`, the Coefficients line (`[v1 ... ] (1x21 double)` for FIR, `[2x6 double]` for SOS), a `Specifications:` block of right-aligned `Name: value` lines and the three `Use fvtool/designfilt/filter` hints. `methods(d)` lists the object functions below.
**Order of work (each step recorded against MATLAB and committed before the next):**
1. The object, its display, `designfilt` for `lowpassfir`/`highpassfir`/`bandpassfir`/`bandstopfir` with `window` (needs `Window`, `CutoffFrequency`, `ScalePassband`, `FilterOrder`; this is `fir1`) and `equiripple`/`ls` with explicit band edges (the `firpm`/`firls` we have), plus the name-value validation and error messages.
2. The IIR responses `lowpassiir`/`highpassiir`/`bandpassiir`/`bandstopiir` with `butter`, `cheby1`, `cheby2` and `ellip` for a given `FilterOrder` (SOS from `zp2sos`, `ScaleSOS`) and the minimum-order forms (`buttord` etc. with `MatchExactly`).
3. The object functions MATLAB lists: `filter`, `filtfilt` (SOS form needed: filtfilt-SOS is still not recovered, so this stays approximate for IIR), `freqz`, `impz`, `impzlength`, `stepz`, `zplane`, `grpdelay`, `phasez`, `phasedelay`, `zerophase`, `filtord`, `firtype`, `info`, `isfir`, `islinphase`, `isminphase`, `ismaxphase`, `isallpass`, `isstable`, `isdouble`, `issingle`, `double`, `single`, `tf`, `zpk`, `ss`, `coeffs`-style property reads, `eq`/`ne`.
4. `lowpass`, `highpass`, `bandpass`, `bandstop` (the signal-processing filter functions), which design with `designfilt` internally (steepness 0.85, 60 dB stopband attenuation, 'auto' impulse response) and so depend on steps 1-3.
Later: `kaiserwin`, `maxflat`, `cls`, `arbmagfir`, `differentiatorfir`, `hilbertfir`, the notch/peak and the DSP System Toolbox responses are not planned (the last need a toolbox that is not installed).
**Open questions to settle by probing MATLAB:** which specification sets `designfilt` accepts per response and method (and its messages for incomplete or inconsistent sets), the exact values it fills in for properties not given, the 'auto' choice between FIR and IIR in `lowpass`/`highpass`, and how the object compares (`eq`).

### designfilt step 1 (done): the object and the four FIR responses
`designfilt` for `lowpassfir`/`highpassfir`/`bandpassfir`/`bandstopfir` builds a `digitalFilter` object (a built-in value class whose property list depends on the specification set; assignments raise MATLAB's "cannot change the specifications" error). Designs: `window` (FilterOrder + CutoffFrequency(1,2); `fir1` with `Window` and `ScalePassband`; an odd-order highpass is the antisymmetric lowpass of cutoff 1-wc modulated by `-(-1)^k`; odd-order bandstop not implemented), `equiripple`/`ls` with FilterOrder + band edges (+ weights; `firpm`/`firls`), and minimum-order `equiripple` from the ripple and attenuation (the smallest order, searched from the `firpmord` estimate, whose design keeps the passband within `dp = (10^(rp/20)-1)/(10^(rp/20)+1)` and the stopband under `10^(-As/20)`). The display is exact except that MATLAB lists the specification properties in an order that changes from one MATLAB session to the next, so MatLivre uses the order the properties were given and the difftest case sorts the display lines. The Coefficients line is shortened with ` ... ] (1xN double)` once the line would pass 179 characters. Not implemented: `kaiserwin`, `maxflat`, `cls`, editing a digitalFilter with `designfilt(d)`, and the other responses.

Odd-order window bandstop (`designfilt('bandstopfir','FilterOrder',31,...,'DesignMethod','window')`): MATLAB returns 32 antisymmetric taps with values near 1e13 (the design is numerically meaningless, with sum 1.84), so it is left unimplemented (MatLivre raises a not-implemented error).

### designfilt step 2 (IIR) - findings so far, not shipped
The denominators are those of `butter`/`cheby1`/`cheby2`/`ellip` in `zp2sos` order 'down' (poles nearest the unit circle last... first section has the larger pole radius), `FilterOrder` for band responses is the total order (6 gives 3 sections). The numerators are the unsolved part: each section is a scaled `[1 2 1]` (lowpass), `[1 -2 1]` (highpass), `[1 0 -1]` (bandpass) or the zeros of the elliptic/chebyshev-II sections, and MATLAB spreads the total gain over the sections in a way that none of `zp2sos`/`ss2sos` with 'up'/'down' and 'none'/'inf'/'two' scaling reproduces. For `butter(4,0.3)` the sections are `[0.1574 0.3148 0.1574]` and `[0.1179 0.2359 0.1179]` (each with unit DC gain); for `cheby1(4,1,0.3)` both have DC gain 0.944 (the square root of the 1 dB ripple); for `ellip(4,1,40,0.3)` the DC gains are 11.8 and 0.0755 (product 0.891). Relative to `ss2sos(...,'down','inf')` the butter sections are scaled by 0.8634 and 0.4303 (product = the leftover gain 0.3717). The rule that produces these scalings is the thing to recover before the IIR designs can match.

### designfilt step 2 (done except the section gains)
`lowpassiir`/`highpassiir`/`bandpassiir`/`bandstopiir` with `butter` (FilterOrder + HalfPowerFrequency), `cheby1` (+ PassbandFrequency, PassbandRipple), `cheby2` (+ StopbandFrequency, StopbandAttenuation), `ellip` (lowpass/highpass: + PassbandFrequency, PassbandRipple, StopbandAttenuation) and the minimum-order forms (`buttord`/`cheb1ord`/`cheb2ord`/`ellipord`, bandpass uses the larger stopband attenuation, bandstop the smaller passband ripple). `FilterOrder` of the band responses is the total order. The object, denominators and the overall transfer function match MATLAB (`sos2tf(d.Coefficients)` on 36 cases, to 1e-9); the **section numerators do not**: MATLAB spreads the total gain over the sections by a rule that fits none of: unit reference gain per section, cumulative L-infinity or L2 scaling of the cascade (with or without the section's own denominator), `zp2sos`/`ss2sos` scalings in either order. MatLivre uses `zp2sos(z,p,k,'down','inf')` with the leftover gain on the last section (3 listed deviations). Fixed-order band `ellip` designs are not accepted by MATLAB with this specification set ("too many parameters ... StopbandAttenuation"), and MatLivre raises the same error.

