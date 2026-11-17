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
| Implemented by this plan (`crates/livre-lang/src/signal/`) | 258 |
| Still to implement and verifiable here (incl. 1 registered without its main syntax) | 8 |
| Not verifiable: requires Wavelet/Audio Toolbox, not installed (see below) | 16 |
| Out of scope (see below) | 34 |

The 258 are checked against MATLAB: see "How it is tested". `isfir`, `tf`, `zpk` and `ss` count as done (digitalFilter methods, which `exist()` cannot see).

Still to do (8): `emd`, `hht`, `fsst`, `kurtogram`, `orderwaveform`, `ordertrack`, `rpmtrack`, and `xwvd` (registered, main syntax not implemented: the cross distribution's rule was not recovered).

## Not verifiable: requires Wavelet/Audio Toolbox, not installed

The MATLAB install used for recording has neither the Wavelet nor the Audio Toolbox, so these cannot be
recorded and are not implemented: `cwt`, `modwt`, `modwpt`, `imodwt`, `wsst`, `wcoherence`, `tqwt`, `cqt`,
`ewt`, `waveletScattering`, `cwtmag2sig`, `wsstridge`, `spectralCrest`, `spectralEntropy`, `spectralFlatness`, `spectralKurtosis`, `spectralSkewness`.
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
- Signal labeling objects and region-of-interest functions (15 functions, object-heavy), and the two feature-extractor frameworks `signalTimeFeatureExtractor` and `signalFrequencyFeatureExtractor` (object-heavy, grouped with the labeling objects).
- EDF file I/O (`edfheader`, `edfinfo`, `edfread`, `edfwrite`): a file-format task.

That removes 34 of the 333 names from the practical target; they stay in `signal.tsv`
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
   outputs), `yulewalk` (all arguments, including `npt` and `lap`).
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
- **`yulewalk`** of high order: the design is ill-conditioned from about order 12, and rounding
  differences in the FFTs and least-squares solves grow to 1e-8 (order 12) and 1e-2 (order 16 with
  a sharp edge). Not recorded, and known to differ: a step at frequency 0 with `lap` 0 or 1
  (rank-deficient equations, MATLAB picks another basic solution) and some complex magnitudes
  (`yulewalk(4,[0 1],[1 1i])`, 6e-4 off in the numerator).
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
- **`filtfilt`** with second-order sections: recovered (see the filtfilt-SOS note below); the IIR design numerators (section gain split, zero-to-pole assignment for elliptic band designs) still differ from MATLAB's.
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
- (`fillgaps` is recovered, see the note at the end.)
- (`filtfilt` with second-order sections is recovered, see the note at the end.)
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

### designfilt step 3 (done): digitalFilter object functions
The class forwards `filter`, `filtfilt`, `fftfilt`, `freqz`, `impz`, `impzlength`, `stepz`, `zplane`, `grpdelay`, `phasez`, `phasedelay`, `zerophase`, `filtord`, `firtype`, `isfir`, `islinphase`, `isminphase`, `ismaxphase`, `isallpass`, `isstable`, `isdouble`, `issingle`, `double`, `single`, `tf`, `zpk`, `ss`, `eq` and `ne` to the numerical functions on the object's numerator/denominator (FIR) or second-order sections (IIR: `filter` is `sosfilt`; `grpdelay`/`phasez`/`phasedelay`/`zerophase` add the section responses, `isminphase`/`ismaxphase` test every section, `zpk` lists the section roots). `filtfilt` of an IIR object raises a not-implemented error (it needs the unrecovered second-order-sections `filtfilt`); `info` is not implemented. 7 listed deviations: FIR `phasedelay` at one point, `ss` of IIR objects (MATLAB's cascade realization), the bandpass group delay at DC, the order of a section pair's zeros, and the 6-digit zeros of a 20th order FIR.

### designfilt step 4 (done, with deviations): lowpass, highpass, bandpass, bandstop
Built on `designfilt(...,'kaiserwin','MinOrder','even')`: cutoff = midpoint of the transition edges, beta from the smallest deviation (Kaiser formula), `'scale'`. Edges: lowpass stop = wp + c*(nyq-wp), highpass stop = wp - c*wp, bandpass likewise on each side, bandstop stops inside the passband gap by c*(gap/2), with c = 0.99 - 0.98*Steepness. The order is the first even n at or above the `kaiserord` estimate whose design keeps the passband peak-to-peak ripple within the ripple and the stopband under the attenuation; FIR is used when the signal length >= 2 x the estimate, with the filter applied by appending N/2 zeros and dropping the first N/2 samples. Deviations: bandpass minimum order is 2 taps higher than MATLAB's in some cases (its stopband check is slightly looser), and the bandstop/IIR branch (SOS plus filtfilt) is not implemented; these cases are in signal_known_deviations.tsv.

### filtfilt with second-order sections (done), and the IIR branch of lowpass/highpass/bandpass/bandstop
`filtfilt(sos,g,x)` (L >= 2 sections; one section is a plain transfer function) filters **section by section**: for each section in turn the signal is padded by odd reflection (2*x(1)-x(i+1)...) with the same pad length for every section, 3 times the order of the whole cascade (`sos2tf` degree, exact trailing zeros dropped: butter(5) gives 15, not 18), started from that section's steady state (`b,a` form), run forward, reversed, started again from the steady state of the reversed start, run again and reversed back; the pad is then cut. The gain is `prod(g(1:min(numel(g),L)))` on the first section (a longer g is an error `There must be at most L+1 scale values.`, an empty g returns `[]`, a g matrix is `Expected input to be a vector.`; non-double arguments give `Invalid data type. Arguments must be double matrices.`; non-finite x gives `Expected input to be finite.`; length <= pad gives `Data length must be larger than N samples.`). This reproduces MATLAB exactly (differences 0 in doubles); it differs from the full-cascade transfer function only at the edges. Ruled out earlier: one filtfilt of the whole cascade with any pad length or initial-state variant. `filtfilt(d,x)` for an IIR digitalFilter is `filtfilt(d.Coefficients,1,x)` (one section: the transfer function of that section).
The IIR branch of the four response functions is `designfilt('<resp>iir', same edges as the FIR design, PassbandRipple 0.1, StopbandAttenuation, 'DesignMethod','ellip')` filtered with `filtfilt`; FIR is used when the signal is at least twice as long as the `kaiserord` estimate. The ellipord bandstop adjustment moves a stopband edge (the other order functions move a passband edge). Design facts recovered: even-order elliptic sections come with the section(s) of the highest-Q poles last (odd order: pure 'down' order, first-order section last); Butterworth/Chebyshev I bandpass sections each take one zero at +1 and one at -1 (`[b 0 -b]`). Not recovered (listed deviations): which zero pair goes with which pole pair inside a pair of band sections of an elliptic design (MATLAB alternates; the edge samples then differ, up to a first sample of -371 for the default bandstop), and the per-section gain split (prefix L-infinity norms were [3.3 34.4 1] and [2.0 2.5 1.1 1], cascade peak exactly 1, last section carries the small gain; no closed form found), which only matters through rounding (4e-9 in an ill-conditioned order-9 case).

### IIR section gains: research queue facts
ellip(6,...) lowpass sections have DC gains 1.57, 14.97 and 0.038; their cumulative L-infinity prefix norms are 3.4, 25.7 and 1; the cascade peak is exactly 1; the first section's gain at DC equals its gain at the passband edge. Ruled out: per-section G^(1/L) at the reference frequency; cumulative L-infinity and L2 scaling on the prefix, on the section's own response and on the internal node 1/A_i; zp2sos/ss2sos up/down with none/inf/two scaling (all miss by 0.04 to 1.2). The gain split only matters through rounding: `filter(d,x)` and `impz(d)` on butter lowpass, butter bandpass and elliptic lowpass objects match MATLAB to 1e-9 as full vectors (cases 8238-8245).

### window, db, scaleFilterSections (done)
`window(@f,N,...)` passes everything to the window function (function handle, or a name as char or string): the output is the first output, two outputs give `Output argument "varargout{2}" (and possibly others) not assigned a value in the execution with "window" function.`; an unknown name gives `Undefined function 'foo' for input arguments of type 'double'.`; anything else `Function to evaluate must be represented as a string scalar, character vector, or function_handle object.` A single-precision N gives a double window (hamming, hann, blackman, rectwin); `hamming`/`hann`/`blackman`/`flattopwin` reject logical and complex N with MATLAB's messages. `db(x,unit,R)`: voltage `10*log10(abs(x).^2./R)` (R default 1; R may be an array or complex, a negative R gives a complex result), power `10*log10(x)` with `Power cannot be negative.` for a negative real; the unit is matched by prefix (`'v'`, `'volt'`, `'p'`, `'POWER'`), in the two-argument form a char is the unit and a number is R; integer and logical input fails in `log10`; single input is computed in double and rounded (MATLAB's single `log10` is 1.3e-6 off, two listed deviations). `scaleFilterSections(N,g)`: g scalar or a vector of L+1 values, section i is multiplied by `g(i)*g(L+1)^(1/L)` (principal root, so a negative g gives complex sections); integer N takes only a scalar g; an N-D array is folded into columns. `hanning(N)` (undocumented legacy): `hann(N+2)` without the end points, `'periodic'` gives `hann(N,'periodic')`, single N gives a single result.

### thd, sfdr, toi, snr, sinad (done for tones; floors are deviations)
All use a power-scaled periodogram of the whole record (nfft = N, Kaiser window beta 38). A tone's power is the sum of the periodogram over its mainlobe divided by the window's equivalent noise bandwidth in bins (`N*sum(w.^2)/sum(w)^2`), and its frequency is the power-weighted centroid of that band; the mainlobe is the monotone-decreasing run of bins either side of the peak, the valley bin on the high side belonging to the next lobe. The fundamental is the largest peak above the DC mainlobe (13 bins; a DC offset is ignored). thd: harmonic k is searched within 2 bins of `k*f1` (aliased folds above Nyquist, otherwise NaN), `r = 10*log10(sum(harmonics 2..n)/fundamental)`, default `n = 6`, frequency scale 1 when `Fs` is omitted; sfdr: the largest local maximum outside the fundamental's lobe and `MSD`; toi: the two largest local maxima are the fundamentals, the products at `2f1-f2` and `2f2-f1`, `OIP3 = mean(Pf + (Pf-Pim)/2)` over the two sides (dB, no 30 dB offset); snr(x,y) is the ratio of summed squared magnitudes; snr(x,Fs,N) divides the fundamental by the periodogram outside the DC guard and the first N harmonic lobes; sinad divides it by everything outside the DC guard and the fundamental. Exact (1e-13) for noise-free tones and for noise when the excluded harmonics are present. Not recovered (listed deviations): harmonics absent from the signal (MATLAB's band for a harmonic that is not there is a single bin at a position that is not `round(k*f1)`, and its power is the bin value itself rather than `sum/enbw`: `hp = P(bin)` for the single-bin case), so the default `N = 6` on a signal with fewer harmonics differs by 0.06 dB in the noise power; floors at -300 dB; tones closer than a mainlobe (toi); single-precision floors. The spectrum-input syntaxes (`'psd'`, `'power'`) are not implemented.

### fillgaps (done)
Two passes over the NaN runs: a forward pass fills the gaps left to right, each from an AR model fitted (Burg, `arburg`) to all samples before the gap with their mean removed, extended recursively over the gap, and the filled values feed the model of the next gap; a backward pass does the same right to left on the reversed samples. Sample `k` of a gap of length `L` is `(1-k/(L+1))*forward + k/(L+1)*backward`; a gap at the start takes only the backward estimate, at the end only the forward one; all-NaN stays NaN; a single non-NaN sample fills everything with it (order 0 = the mean). `MAXLEN` keeps the last (first) `MAXLEN` samples before (after) the gap and an `ORDER` larger than `MAXLEN` is an error. The automatic order minimises `N*log(e) + 2*q` over `q = 1..N-1` (`e` = the Burg error variance, `N` the number of samples used, picked separately for the two directions); a fixed or infinite order is truncated to `N-1`. Single, integer, logical, N-D input and Inf are errors (`Expected X to be one of these types: double`). Exact (1e-13) except: two gaps with automatic order when the backward data contain the AR extension of a later gap (the AIC picks the maximal order, a degenerate fit), and complex input.

### stftmag2sig (done for 'gla')
Griffin-Lim exactly as MATLAB (1e-13): start `Sest = S.*exp(i*phase)` (zero phase), then each iteration is `x = istft(Sest, ConjugateSymmetric for centered/twosided); S2 = stft(x); Sest = S.*exp(i*angle(S2))`; it stops after `MaxIterations` (default 100) or at the first iteration whose spectral convergence `norm(abs(S2)-S)/norm(S)` is below `InconsistencyTolerance` (default 1e-4; the documented inconsistency `norm(S2-Sest)/norm(Sest)` does not reproduce the stopping point); the result is `istft(Sest)` after the last iteration (`MaxIterations=0` is `istft(S)`). Window default `hann(128,'periodic')`, overlap `floor(0.75*numel(window))`, FFT length at least the window length, `FrequencyRange` centered/onesided/twosided, 3-D input is a multichannel signal, `InputTimeDimension='downrows'` transposes. Outputs: `t = (0:length(x)-1)'/Fs` (or `*Ts` as a duration, samples when omitted) and `info` with `ExitFlag` (1 when the tolerance was met), `NumIterations`, `Inconsistency` (the last spectral convergence, `Inf` for zero iterations), `ReconstructedSTFT = Sest` and `ReconstructedPhase = angle(Sest)`. `Method='gd'` fails with MATLAB's `Undefined function 'dlarray'` (Deep Learning Toolbox). Not recovered: `'fgla'` (Perraudin's momentum step matches MATLAB after one iteration only), `'legla'`, `InitializePhaseMethod='random'`, a nonzero `InitialPhase` (a constant offset has no effect in MATLAB), and MATLAB's single-precision iterations.

### strips (done)
One line per column made of NaN-terminated strips of N samples (the last strip padded with NaN; the NaNs are negative NaN), the first strip on top: `y = x*scale - mid + offset`, `mid = (max+min)/2` over all columns, strip spacing `d = 1.25*range` (1 for a constant signal, then `YLim = [0 S*d]`), `YLim = [-0.75*range, S*d - 0.5*d + 0.1*d]` (the operation order matters to the last bit), `YTick` = the offsets, `YTickLabel = num2str` of the strip start samples (or seconds for `strips(x,sd,fs)`, where N = round(sd*fs)), `XLim = [0 N]` or `[0 sd]`; a complex X warns `X vector complex.  I'm ignoring the imaginary part.`. The recorded cases use rounded data (sin/cos differ in the last bit).

Absent-harmonic data point (thd/snr): with f1 = 77.3 Hz, N = 1024, fs = 1000, MATLAB's single-bin lobes are at 0-based bins 395 (k=5, expected 395.84) and 476 (k=6, expected 474.93; bin 477 is the local maximum); `round`, `floor` and a biased centroid estimate of f1 each fit at most one of them, and the lobe power is the raw bin value.

### kurtogram (closed: not recovered, two rounds)
Known: the call is `kurtogram(x, Fs, level)` (the second argument is the sample rate, so `kurtogram(x,[],0)`); the level-0 row is exactly the centred real kurtosis minus 3, `mean((x-mean(x)).^4)/mean((x-mean(x)).^2)^2 - 3`, constant across the 3 columns (1.94393545547285 for `x = sin(0.9n)+0.4 sin(0.05 n^1.3)+0.2 cos(2.7n)+3*(mod(n,37)==0)`, n = 1:256); for `level = 1` the matrix is 2 x 6 and the second row is `[K_low K_low K_low K_high K_high K_high]` = `[0.200825242192111 x3, 5.03012351574369 x3]` for that signal, so each band's value is repeated over its share of the `3*2^level` columns; the default level is `max(0, floor(log2(N)-6))`; W holds the equivalent window sizes `2^(n+1)`. Ruled out for the level-1 values: Antoni's prototype `h = fir1(16,0.4).*exp(2i*pi*(0:16)*0.125)` with `g = h(1+mod(1-n,16)).*(-1).^(1-n)` (n = 2:16) applied as filter/conv, decimating from sample 1 or 2, with or without trailing zeros or `same` trimming, and the complex kurtosis `mean|z|^4/mean|z|^2^2 - 2` (best 0.2084 / 5.5456, 0.0218 / 6.668, -0.0346 / 6.463). The exact padding, decimation phase, filter lengths per level and any kurtosis correction remain unknown.

fsst, second round: re-implemented the rule above for a 96-sample signal with `kaiser(16,10)` (output 9 x 96): every cell is within 1e-3 relative of MATLAB's magnitude after a common scale of 1.0008, so the squeezing positions agree and the remaining difference is a smooth amplitude one (window scaling or segment alignment), not the rounding-boundary cells; not pursued further.

fsst, correction: the common scale of 1.0008 was an artefact of a least-squares fit. Cell by cell the ratio mine/MATLAB is exactly 1 (1e-15) for most cells, in particular rows 7-9 and most columns, and differs only in rows 1-6 of some columns (up to 8 % and 20 % in rows 4 and 5, residual 1.6 % after the best scale), none of the window sums (symmetric 5.8670, periodic 6.2578, w(L/2+1) 0.9791, norm 2.0577) gives it. So the differences are individual cells whose squeezed position or contributions differ from MATLAB's (the earlier rounding-boundary diagnosis stands); `fsst` stays closed as not recovered.

### Vold-Kalman family (orderwaveform, ordertrack, rpmtrack): closed, blocked by the solver
For constant rpm the first-order filter is `|H(df)| = 1/(1+(df/a)^4)` with `a` about 0.633*Bandwidth (bw 2: a = 1.266; bw 5: 3.166; bw 10: 6.33; the exact half-power value would be 0.623*bw), equivalent to the second-difference weight `r = (fs/(2*pi*a))^4`. MATLAB solves it with `pcg` (`signal.internal.vk` warns that the tolerance may not be achievable and that PCG did not converge for FilterOrder 2 and narrow bandwidths) and the output at df = 0 is 1.00003 / 1.0006 / 1.00997 for bw 2 / 5 / 10 instead of 1, so its results carry an iteration error of 1e-3 to 1e-2 that a direct solver cannot reproduce to 1e-9.

Vold-Kalman, second check: the model is the system `(I + r*D'*D) x = 2*conj(c).*y` with `D = diff(speye(N),2)`, carrier `c = exp(2i*pi*f0*t)`, output `real(x.*c)`; a direct solve matches `orderwaveform` (constant rpm, FilterOrder 1, Bandwidth 5, N = 3000) to 1.2e-5 with `a/bw = 0.632911` (`r = (fs/(2*pi*a))^4`; 4.6e-4 with 0.6332), a constant not explained by the half-power formula (0.6232). Plain `pcg(A,b,tol,maxit)` without a preconditioner returns the zero start (flag 1, 0 iterations) for this r and N, so MATLAB's solver settings are different and not inferable.
