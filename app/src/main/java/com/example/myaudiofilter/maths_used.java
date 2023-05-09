package com.example.myaudiofilter;

import android.content.Context;
import android.widget.Toast;

import java.util.Arrays;


public class maths_used {
    static double[] values1;

    public static double sinc(double x) {
        if (x == 0) {
            return 1.0;
        } else {
            return Math.sin(Math.PI * x) / (Math.PI * x);
        }
    }

    public static double window(int i, int N) {
        return 0.42 - 0.5 * Math.cos(2 * Math.PI * i / (N - 1)) + 0.08 * Math.cos(4 * Math.PI * i / (N - 1));
    }

    public double getWindowValue(int n, int filterLength, int windowType) {
        double value;
        switch (windowType) {
            case 0: // Rectangular window
                value = 1.0;
                break;
            case 1: // Bartlett (triangular) window
                value = 1.0 - 2.0 * Math.abs(n - (filterLength - 1) / 2.0) / (filterLength - 1);
                break;
            case 2: // Hann (Hanning) window
                value = 0.5 - 0.5 * Math.cos(2 * Math.PI * n / (filterLength - 1));
                break;
            case 3: // Hamming window
                value = 0.54 - 0.46 * Math.cos(2 * Math.PI * n / (filterLength - 1));
                break;
            default:
                throw new IllegalArgumentException("Invalid window type: " + windowType);
        }
        return value;
    }

    public static double[] convolve(int[] x, double[] h) {

        // Compute length of output array
        int yLength = x.length + h.length - 1;

        // Create output array
        double[] y = new double[yLength];

        // Compute convolution
        for (int i = 0; i < yLength; i++) {
            for (int j = 0; j < h.length; j++) {
                if (i - j >= 0 && i - j < x.length) {
                    y[i] += x[i - j] * h[j];
                }
            }
        }
        return y;

    }

    private static byte[] toByteArray(int[] data) {
        byte[] bytes = new byte[data.length * 2];
        for (int i = 0; i < data.length; i++) {
            bytes[i * 2] = (byte) (data[i] & 0xff);
            bytes[i * 2 + 1] = (byte) ((data[i] >> 8) & 0xff);
        }
        return bytes;
    }

    public static int[] convertByteArrayToIntArray(byte[] audioData, int bytesPerSample) {
        int[] samples = new int[audioData.length / bytesPerSample];
        int sampleIndex = 0;
        for (int i = 0; i < audioData.length; i += bytesPerSample) {
            int sample = 0;
            for (int j = 0; j < bytesPerSample; j++) {
                sample |= (audioData[i + j] & 0xFF) << (j * 8);
            }
            samples[sampleIndex++] = (short) sample;
        }
        return samples;
    }

    public static byte[] inttoBytes(double[] filteredSignal, int bytesPerSample) {
        byte[] newData = new byte[filteredSignal.length * bytesPerSample];
        int dataIndex = 0;
        for (int i = 0; i < filteredSignal.length; i++) {
            int sample = (int) filteredSignal[i];
            for (int j = 0; j < bytesPerSample; j++) {
                newData[dataIndex++] = (byte) (sample & 0xFF);
                sample >>= 8;
            }
        }
        return newData;
    }

    public static byte[] interleaveData(byte[] leftChannelData, byte[] rightChannelData) {
        byte[] interleavedData = new byte[leftChannelData.length + rightChannelData.length];
        int i = 0, j = 0, k = 0;
        while (i < leftChannelData.length && j < rightChannelData.length) {
            interleavedData[k++] = leftChannelData[i++];
            interleavedData[k++] = leftChannelData[i++];
            interleavedData[k++] = rightChannelData[j++];
            interleavedData[k++] = rightChannelData[j++];
        }
        while (i < leftChannelData.length) {
            interleavedData[k++] = leftChannelData[i++];
            interleavedData[k++] = leftChannelData[i++];
        }
        while (j < rightChannelData.length) {
            interleavedData[k++] = rightChannelData[j++];
            interleavedData[k++] = rightChannelData[j++];
        }
        return interleavedData;
    }

    public static double[] computeFilterTaps(int numTaps, int numBands, double[] bands,
                                             double[] des, double[] weight, Remez.FilterType type) {

        // Set up input parameters
        final double[] x = new double[numTaps];
        final double[] y = new double[numTaps];
        final double[] ad = new double[numTaps];
        final double[] h = new double[numTaps];
        //System.out.println(" initial h"+Arrays.toString(h));
        // Call Remez function
        boolean success = Remez.remez(numTaps, numBands, bands, des, weight, type, x, y, ad, h);

        if (!success) {
            System.err.println("Error: Failed to converge");
            return null;
        }
        System.out.println(" final h" + Arrays.toString(h));
        System.out.println("Sucess");
//        // Get filter taps
//        int numTapsHalf = (numTaps + 1) / 2;
//        double[] taps = new double[numTapsHalf];
//        System.arraycopy(h, 0, taps, 0, numTapsHalf);
//        for (int i = 0; i < numTapsHalf; i++) {
//            taps[i] *= 2.0;
//        }
//        if (numTaps % 2 == 0) {
//            taps[numTapsHalf - 1] /= 2.0;
//        }

        return h;
    }

    public static double[] computeRemezBands(String RemezBandsSring) {
        float value1, value2, value3, value4, value5, value6;
        String[] values = RemezBandsSring.split(",");

        try {

            //arr[5] = 10; // This will throw an ArrayIndexOutOfBoundsException

            if (values.length == 4) {
                value1 = Float.parseFloat(values[0]);
                value2 = Float.parseFloat(values[1]);
                value3 = Float.parseFloat(values[2]);
                value4 = Float.parseFloat(values[3]);
                value5 = 0.0f; // Default value for value5
                value6 = 0.0f; // Default value for value6
            }
            else if  (values.length == 6) {
                value1 = Float.parseFloat(values[0]);
                value2 = Float.parseFloat(values[1]);
                value3 = Float.parseFloat(values[2]);
                value4 = Float.parseFloat(values[3]);
                value5 = Float.parseFloat(values[4]);
                ; // Default value for value5
                value6 = Float.parseFloat(values[5]);
            } else {
                value1 = 0.0f; // Default value for value1
                value2 = 0.2f; // Default value for value2
                value3 = 0.25f; // Default value for value3
                value4 = 1f; // Default value for value4
                value5 = 0.0f; // Default value for value5
                value6 = 0.0f; // Default value for value6
            }
            double[] values1 = new double[]{value1, value2, value3, value4, value5, value6};
            if (values.length == 4) {
                return new double[]{value1, value2, value3, value4};
            } else {
                return values1;
            }

        } catch (Exception e) {
            values1 = new double[]{0.0, 0.2, 0.25, 1};
            return values1;
        }

    }

    public static double[] computeRemezWeights(String RemezBandsSring) {
        int value1, value2, value3, value4, value5, value6;
        String[] values = RemezBandsSring.split(",");

        try {
            if (values.length == 2) {
                value1 = Integer.parseInt(values[0]);
                value2 = Integer.parseInt(values[1]);
//                value3 = Integer.parseInt(values[2]);
//                value4 = Integer.parseInt(values[3]);
//                if (values.length >= 5) {
//                    value5 = Integer.parseInt(values[4]);
//                } else {
//                    value5 = 0; // Default value for value5
//                }
                if (values.length == 3) {
                    value6 = Integer.parseInt(values[2]);
                } else {
                    value6 = 0; // Default value for value6
                }
            } else {
                value1 = 1; // Default value for value1
                value2 = 0; // Default value for value2
                value6 = 0; // Default value for value3
//                value4 = 0; // Default value for value4
//                value5 = 0; // Default value for value5
//                value6 = 0; // Default value for value6
                //value3, value4, value5,
            }
            double[] values1 = new double[]{value1, value2, value6};
            if (values.length == 2) {
                return new double[]{value1, value2};
            } else {
                return values1;
            }

        } catch (Exception e) {
            values1 = new double[]{0,1,0};
            return values1;
        }

    }
}

