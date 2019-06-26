"""This class will simulate a calcium dataset.

Author: s w keemink
"""
from __future__ import division
import numpy as np
from . import simcalclib as sclib
import tifffile


class SimCalc():
    """Simlates a calcium imaging dataset given parameters.

    Parameters
    ----------
    stimprop : float
        proportion of cells that are tuned to the stimulus
    h : int
        Height and width of image in pixels
    N : int
        Number of neurons
    T : float
        Stimulus time. Total time will be StimCycles*2*T
    StimCycles : int, optional (default: 4)
        How many stimulus cycles within T
    dt : float, optional (default: 0.01)
        Simulation time
    """

    def __init__(self, h, N, T, StimCycles=4, dt=0.01):
        """Initialisation."""
        self.h = h
        self.N = N
        self.T = T
        self.StimCycles = StimCycles
        self.dt = dt

    def gen_calcium(self, Tau_r, Tau_d, A, p, rates=None, sp=None):
        """Generation of the spikes.

        Spikes are generated according to Poisson statistics.

        Calcium is first modelled with a rise and fall time, then passed
        through a polynomial to model nonlinearities:
        c_d' = -c_d/Tau_d + s(t)
        c_r' = -c_r/Tau_r + s(t)
        c = c_d + c_r

        F(c) = 1 + A[c+p2(c^2-c)) + p3(c^3-c)]

        Parameters
        ----------
        Tau_r, Tau_d : floats
            Rise and decay times of calcium
        A : float/array
            Spike response magnitude.
            If Array, should be of shape (N,1), to give every neuron a
            different value.
        p : array/list
            Calcium dynamics polynomial variables (p = [p2,p3])
        rates : array, optional (default: None)
            If None, the firing rates follow a normal distribution with
            a given mean and sigma. If given, should be an array giving the
            firing rates for every neuron.
        sp : array, optional (default: None)
            If None, will make spikes according to parameters. Otherwise,
            will use provided spikes to make calcium trace.
        """
        N = self.N
        dt = self.dt
        if sp is None:
            # get basic variables
            T = self.T
            Cycles = self.StimCycles
            nFrames = int(T / dt) * 2 * Cycles
            times = np.arange(0, T * 2 * Cycles, dt)

            # set up firing rates
            if rates is None:
                rates = np.random.normal(0.5, 1, N)
                rates[rates < 0] = 0

            # Model Poisson firing
            sp = np.zeros((N, nFrames))
            for i in range(N):
                sp[i, :] = sclib.spike_data_poisson(
                    T, dt, [rates[i], rates[i] * 2] * Cycles, trials=1)
        else:
            nFrames = sp.shape[1]

        # Simulate Calcium dynamics
        c = np.zeros((N, nFrames))  # calcium traces
        F = np.zeros(c.shape)
        c_d = np.zeros(N)  # calcium decay dynamics
        c_r = np.zeros(N)  # calcium rise dynamics
        for i, t in enumerate(times):
            # iterate calcium levels
            c_d += dt * (sp[:, i] / dt - c_d / Tau_d)
            c_r += dt * (sp[:, i] / dt - c_r / Tau_r)
            c[:, i] = c_d - c_r

        # calculate basic dye levels
        # F = A*c

        maxc = (-2 * p[0] - np.sqrt(4 * p[0]**2 +
                                    12 * p[1] * (sum(p) - 1))) / (6 * p[1])
        maxf = (A * (maxc + p[0] * (maxc**2 - maxc) + p[1] * (maxc**3 - maxc)))
        F = A * (c + p[0] * (c**2 - c) + p[1] * (c**3 - c))
        if np.sum(c > maxc) > 0:
            F[c > maxc] = maxf

        self.sp = sp
        self.c = c
        self.F = F

    def gen_spat_kernels(self, locs=None, sizes=None, covs=None,
                         use_rings=None):
        """Generate the spatial kernels.

        Parameters
        ----------
        locs : array (number of neurons by 2), optional
            Locations of every neuron. If None, just randomizes.
        sizes : array, optional
            Sizes of every neuron. If None, just randomizes.
        covs : array, optional
            The covariance between width and height
            of each neuron. If None, sets all to 0.
        use_rings : array of bools, optional
            Whether to use Ring neurons or not. All True by default.

        """
        # get basic variables
        N = self.N
        h = self.h
        spat_kernels = [0] * N
        masks = [0] * N

        # get locations and sizes
        if locs is None:
            locs = np.random.randint(0, h, (N, 2))
        if sizes is None:
            sizes = np.random.normal(h / 2, np.sqrt(h / 2), size=N)
        if covs is None:
            covs = np.zeros(N)
        if use_rings is None or use_rings:
            use_rings = np.ones(N, dtype=bool)

        # generate shapes
        for i in range(N):
            cov = np.array([[sizes[i], covs[i]], [covs[i], sizes[i]]])
            # make spatial kernel
            spat_kernels[i] = sclib.makeGaussianFilter(
                locs[i, 0], locs[i, 1], cov, h, ring=use_rings[i])
            spat_kernels[i] /= spat_kernels[i].max()

            # make ROI mask
            masks[i] = spat_kernels[i] > 0.5

            # add offset to original kernel
            spat_kernels[i] += masks[i] / 5

        self.spat_kernels = spat_kernels
        self.masks = masks

    def gen_npil(self, eta=0.05, F0=1, nFilts=10):
        """Generation of neuropil, both its mask and responses.

        Temporal component:
        B' = eta*dW

        Spatial component:
        Sum of 2d-Gaussians with random locations and widths.
        The spatial component is normalized so the mean is 1.

        Parameters
        ----------
        eta : float
            Magnitude and timescale of brownian motion.
        F0 : float
            The drifting baseline
        nFilts : int
            Number of background filters to use for spatial component.
            Set to 0 to have flat background.
        """
        # get basic variables
        h = self.h
        T = self.T
        dt = self.dt
        Cycles = self.StimCycles
        nFrames = int(T / dt) * 2 * Cycles

        # run brownian motion for temporal component
        drive = np.concatenate([np.ones(int(T / dt)),
                                np.ones(int(T / dt)) * 1.05] * Cycles) * F0
        Bg = drive + eta * \
            np.cumsum(np.random.normal(size=nFrames)) * np.sqrt(dt)
        Bg[Bg < 0] = 0

        # start spatial component
        Bg_img = np.ones((1, h, h))

        if nFilts > 0:
            # random gaussian locations
            x, y = (np.random.randint(0, h, nFilts),
                    np.random.randint(0, h, nFilts))
            # loop over locations
            for i in range(len(x)):
                cov = [[np.random.normal(h * 2, np.sqrt(h * 2)), 0],
                       [0, np.random.normal(h * 2, np.sqrt(h * 2))]]

                mask_temp = sclib.makeGaussianFilter(x[i], y[i], cov, h)
                Bg_img[0, :, :] += mask_temp

            # normalize so mean is 1
            Bg_img /= Bg_img.mean()

        # make video
        video = np.repeat(Bg_img, nFrames, axis=0) * Bg[:, None, None]

        self.Bg_trace = Bg
        self.Bg_img = Bg_img[0, :, :]
        self.npil = video

    def gen_data(self, p=None, ampl=1):
        """Generation of final data, combining spikes, masks and neuropil.

        Run gen_npil, gen_spat_kernels, and gen_calcium first.

        Parameters
        ----------
        p : array/list
            Calcium dynamics polynomial variables (p = [p2,p3])
            not used at the moment
        ampl : float, optional
            Multiplication of output
        """
        # get basic variables
        N = self.N
        h = self.h
        T = self.T
        dt = self.dt
        Cycles = self.StimCycles
        nFrames = int(T / dt) * 2 * Cycles
        F = self.F
        kernels = self.spat_kernels

        # start cell video
        data = np.ones((nFrames, h, h))
        print('Of ' + str(N) + ' cells:')
        for cell in range(N):
            print(cell)
            cell_img = kernels[cell].reshape((1, h, h))
            data += np.repeat(cell_img, nFrames, axis=0) * \
                F[cell, :, None, None]
            # TODO: speed up above two lines

        # add neuropil
        data += self.npil - 1
        # data -= 1
        # data = (self.npil)*(1+data)-1
        # data[data<0]=0
        # data = data*ampl

        # optionally apply nonlinearity
        # data = 1 + (data + p[0]*(data**2 - data) + p[1]*(data**3 - data))
        data[data < 0] = 0

        # measure video through poisson
        video = np.random.poisson(data)

        self.data = data
        self.video = video

    def show_video(self, roi_id, downt):
        """Showing Holoviews video of data.

        Parameters
        ----------
        roi_id : int
            Which ROI to show
        downt : int
            Downsampling in time

        Returns
        -------
        HoloMap
            Video
        HoloMap
            Frame indicator to put on top of Curves
        """
        # load tiff
        data = self.video[:, ::-1, :]

        # make downsampled video
        numFrames = data.shape[0]
        bnds = (0, 0, data.shape[1], data.shape[2])
        frames = {f: hv.Image(np.mean(data[f:f + downt, :, :], axis=0),
                              bounds=bnds) for f in range(0, numFrames, downt)}
        vid = hv.HoloMap(frames, kdims=['frames'])

        # make frame indicator
        vlines = {f: hv.VLine(f + downt / 2)
                  for f in range(0, numFrames, downt)}
        frame_indicator = hv.HoloMap(vlines, kdims=['frames'])

        return vid, frame_indicator

    def gen_masks(self, threshold=0.5):
        """Generate the masks for each ROI.

        Parameters
        ----------
        threshold : float, optional
            Which threshold to use for mask generation.

        Returns
        -------
        list
            list of binary arrays representing the masks for each ROI
        """
        # get basic variables
        N = self.N
        masks = [0] * N

        # turn cell kernels into masks
        for i in range(N):
            # make ROI mask
            masks[i] = self.spat_kernels[i] > threshold

        return masks

    def save_tiff(self, fname='data.tif', dtype=np.int16):
        """Saving the generated data as a tiff file."""
        # shpe = self.video.shape
        # video = np.zeros((shpe[2], self.h, self.h), dtype=dtype)
        # for i in range(shpe[2]):
        #     video[i, :, :] = self.video.astype(dtype)[:, :, i]
        video = self.video.astype(dtype)
        tifffile.imsave(fname, video, imagej=True)
