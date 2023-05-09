package com.example.myaudiofilter;

import static com.example.myaudiofilter.maths_used.computeFilterTaps;
import static com.example.myaudiofilter.maths_used.computeRemezBands;

import android.annotation.SuppressLint;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.graphics.Color;
import android.net.Uri;
import android.os.Bundle;
import android.Manifest;


import androidx.activity.result.ActivityResultLauncher;
import androidx.activity.result.contract.ActivityResultContracts;
import androidx.annotation.NonNull;
import androidx.annotation.Nullable;
import androidx.appcompat.app.AppCompatActivity;

import android.provider.DocumentsContract;
import android.view.View;

import androidx.core.app.ActivityCompat;
import androidx.core.content.ContextCompat;
import androidx.core.content.FileProvider;


import com.example.myaudiofilter.databinding.ActivityMainBinding;

import android.view.Menu;
import android.view.MenuItem;
import android.widget.AdapterView;
import android.widget.Button;
import android.widget.EditText;
import android.widget.Spinner;
import android.widget.TextView;
import android.widget.Toast;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import com.example.myaudiofilter.maths_used;

public class MainActivity extends AppCompatActivity {

    //private static final int REQUEST_CODE_LOAD_FILE = 1;
    private static final int MY_PERMISSIONS_REQUEST_WRITE_EXTERNAL_STORAGE = 1;
    static private byte[] audioData;
    static private byte[] newData;
    static byte[] header = new byte[44];
    static int[] audioData_Decimal;
    static  double[] filterCoefficients;
    int bytesPerSample;
    int selectedWindow;
    static int numChannels ;
    static int bitsPerSample ;
    static int[] samples;
    static double[] filteredSignal;
    static double[] bands;
    static double[] des ;
    static double[] weight;
    static int numBands;
    static int numTaps;
    int[] leftChannelSamples;
    int[] rightChannelSamples;
    static double[] filteredLeftChannel;
    static double[] filteredRightChannel;
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
// Check for the WRITE_EXTERNAL_STORAGE permission
        if (ContextCompat.checkSelfPermission(this, Manifest.permission.WRITE_EXTERNAL_STORAGE) != PackageManager.PERMISSION_GRANTED) {
            // Request the permission if it has not been granted
            ActivityCompat.requestPermissions(this, new String[]{Manifest.permission.WRITE_EXTERNAL_STORAGE}, MY_PERMISSIONS_REQUEST_WRITE_EXTERNAL_STORAGE);
        }

        Button loadButton = findViewById(R.id.loadFile);
        Button applyFilterButton = findViewById(R.id.ReduceAmplitude);
        Button saveButton = findViewById(R.id.saveFile);

        TextView BitsPerSample = findViewById(R.id.BitsperSample);
        EditText RemezBands = findViewById(R.id.remezCoff);
        EditText RemezWeights = findViewById(R.id.remezCoffWeights);
        EditText RemezTaps = findViewById(R.id.NumTaps);
        TextView NumChannels = findViewById(R.id.NumChannel);
        TextView selectWindow = findViewById(R.id.selectWindow);
        Spinner filterSpinner = findViewById(R.id.filterspinner);
        filterSpinner.setSelection(-1);

        loadButton.setOnClickListener(view -> {
            Intent intent = new Intent(Intent.ACTION_OPEN_DOCUMENT);
           // intent.putExtra(DocumentsContract.EXTRA_INITIAL_URI, uriToLoad);
            intent.addCategory(Intent.CATEGORY_OPENABLE);
            intent.setType("*/*");
            someActivityResultLauncher.launch(intent);
        });
        filterSpinner.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {
            public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {
                Toast.makeText(MainActivity.this, filterSpinner.getSelectedItem().toString() + "Selected", Toast.LENGTH_SHORT).show();
                String selectedValue = filterSpinner.getSelectedItem().toString();
                //selectedWindow = Integer.parseInt(selectedValue);
                selectedWindow=position;
            }
            @Override
            public void onNothingSelected(AdapterView<?> adapterView) {
                Toast.makeText(MainActivity.this, "Nothing is Selected", Toast.LENGTH_SHORT).show();
            }
        });
        applyFilterButton.setOnClickListener(view -> {
            //System.out.println("Input Data" + audioData);
            // Convert the input signal to an array of 16 bit integer
            if (numChannels == 1 && bitsPerSample == 8) {
                // mono, 8-bit audio data
                System.out.println("Still have to do processing");

            } else if (numChannels == 1 && bitsPerSample == 16) {
                // mono, 16-bit audio data
                samples = maths_used.convertByteArrayToIntArray(audioData, bytesPerSample);
                System.out.println(("Samples Length :" + samples.length));

                // samples = maths_used.convertByteArrayToIntArray(audioData, bytesPerSample);

// Define the parameters of the low-pass filter
                double samplingRate = 44100.0; // Hz
                double cutoffFrequency = 2000.0; // Hz
                double filterLengthSeconds = 0.01; // seconds
                maths_used mathsUsed = new maths_used();
// Compute the filter coefficients using a suitable method

                if (selectedWindow == 4) {
                    String RemezBandsString = RemezBands.getText().toString();
                    String RemezWeightsString = RemezWeights.getText().toString();
                    try {
                        String NumberRemezTaps = RemezTaps.getText().toString();
                        numTaps = Integer.parseInt(NumberRemezTaps);
                    } catch (Exception e) {
                      numTaps= 19;
           }
                    double [] values=maths_used.computeRemezBands(RemezBandsString);
                    double [] Weight_values=maths_used.computeRemezWeights(RemezWeightsString);
                    System.out.println(Weight_values.length);

                    //numBands = 2;

                    Remez.FilterType type = Remez.FilterType.Bandpass;
//                    double Fc = 1000.0;
//                    double fs = 8000.0;
////                    bands[1] = 2.0 * Fc / fs;
////                    bands[2] = bands[1];

                        if (values.length == 4) {
                            numBands = 2;
                            weight = new double[]{1.0, 1.0};
                            bands = new double[]{values[0], values[1], values[2], values[3]};
                            des= new double[]{Weight_values[0], Weight_values[1]};
                        }
                        else if (values.length == 6){
                            numBands = 3;
                            weight = new double[]{1.0, 1.0, 1.0};
                            bands = new double[]{values[0],values[1], values[2], values[3], values[4], values[5]};
                            des= new double[]{Weight_values[0], Weight_values[1], Weight_values[2]};
                        }
                    System.out.println("Weight"+Arrays.toString(weight));
                    System.out.println("Bands"+Arrays.toString(bands));
                    System.out.println("Desired Values"+Arrays.toString(des));
                   try {
                        filterCoefficients = computeFilterTaps(numTaps, numBands, bands, des, weight, type);
                    } catch (Exception e) {
                        Toast.makeText(MainActivity.this, "Change Remez Bands between 0,1", Toast.LENGTH_SHORT).show();
                    }
                    //
                   // filterCoefficients = computeFilterTaps(numTaps, numBands, bands, des, weight, type);
                } else {
                    int filterLengthSamples = (int) Math.ceil(filterLengthSeconds * samplingRate);
                    filterCoefficients = new double[filterLengthSamples];
                    double alpha = 2 * Math.PI * cutoffFrequency / samplingRate;
                    for (int n = 0; n < filterLengthSamples; n++) {
                        double hn = alpha * maths_used.sinc(alpha * (n - (filterLengthSamples - 1) / 2.0));
                        filterCoefficients[n] = hn * mathsUsed.getWindowValue(n, filterLengthSamples, selectedWindow);
                    }
                    double sum = 0.0;
                    for (int i = 0; i < filterLengthSamples; i++) {
                        sum += filterCoefficients[i];
                    }
                    for (int i = 0; i < filterLengthSamples; i++) {
                        filterCoefficients[i] /= sum;
                    }
                }

                try {
                    filteredSignal = maths_used.convolve(samples, filterCoefficients);// Apply the filter to the input signal using convolution
                    newData = maths_used.inttoBytes(filteredSignal, bytesPerSample);// Convert the filtered signal back to bytes
                    Toast.makeText(MainActivity.this, "Filter Applied", Toast.LENGTH_SHORT).show();
                } catch (Exception e) {
                    Toast.makeText(MainActivity.this, "Try Changing Remez Bands ", Toast.LENGTH_SHORT).show();
                }
            } else if(numChannels == 2 && bitsPerSample == 16) {
                byte[] leftChannel = new byte[audioData.length / 2];
                byte[] rightChannel = new byte[audioData.length / 2];
                for (int i = 0; i < audioData.length; i += bytesPerSample * numChannels) {
                    for (int j = 0; j < bytesPerSample; j++) {
                        leftChannel[(i / (bytesPerSample * numChannels)) * bytesPerSample + j] = audioData[i + j];
                        rightChannel[(i / (bytesPerSample * numChannels)) * bytesPerSample + j] = audioData[i + j + bytesPerSample];
                    }
                }
                leftChannelSamples = maths_used.convertByteArrayToIntArray(leftChannel,bytesPerSample);
                rightChannelSamples = maths_used.convertByteArrayToIntArray(rightChannel,bytesPerSample);
                double samplingRate = 44100.0; // Hz
                double cutoffFrequency = 1000.0; // Hz
                double filterLengthSeconds = 0.01; // seconds
                maths_used mathsUsed = new maths_used();
                // stereo, 16-bit audio data
                if (selectedWindow == 4) {
                    String RemezBandsString = RemezBands.getText().toString();
                    String RemezWeightsString = RemezWeights.getText().toString();
                    try {
                        String NumberRemezTaps = RemezTaps.getText().toString();
                        numTaps = Integer.parseInt(NumberRemezTaps);
                    } catch (Exception e) {
                        numTaps= 19;
                    }
                    double [] values=maths_used.computeRemezBands(RemezBandsString);
                    double [] Weight_values=maths_used.computeRemezWeights(RemezWeightsString);
                    System.out.println(Weight_values.length);

                    //numBands = 2;

                    Remez.FilterType type = Remez.FilterType.Bandpass;
//                    double Fc = 1000.0;
//                    double fs = 8000.0;
////                    bands[1] = 2.0 * Fc / fs;
////                    bands[2] = bands[1];

                    if (values.length == 4) {
                        numBands = 2;
                        weight = new double[]{1.0, 1.0};
                        bands = new double[]{values[0], values[1], values[2], values[3]};
                        des= new double[]{Weight_values[0], Weight_values[1]};
                    }
                    else if (values.length == 6){
                        numBands = 3;
                        weight = new double[]{1.0, 1.0, 1.0};
                        bands = new double[]{values[0],values[1], values[2], values[3], values[4], values[5]};
                        des= new double[]{Weight_values[0], Weight_values[1], Weight_values[2]};
                    }
                    System.out.println("Weight"+Arrays.toString(weight));
                    System.out.println("Bands"+Arrays.toString(bands));
                    System.out.println("Desired Values"+Arrays.toString(des));
                    try {
                        filterCoefficients = computeFilterTaps(numTaps, numBands, bands, des, weight, type);
                    } catch (Exception e) {
                        Toast.makeText(MainActivity.this, "Change Remez Bands between 0,1", Toast.LENGTH_SHORT).show();
                    }
                   // filterCoefficients = computeFilterTaps(numTaps, numBands, bands, des, weight, type);
                } else {
                    int filterLengthSamples = (int) Math.ceil(filterLengthSeconds * samplingRate);
                    filterCoefficients = new double[filterLengthSamples];
                    double alpha = 2 * Math.PI * cutoffFrequency / samplingRate;
                    for (int n = 0; n < filterLengthSamples; n++) {
                        double hn = alpha * maths_used.sinc(alpha * (n - (filterLengthSamples - 1) / 2.0));
                        filterCoefficients[n] = hn * mathsUsed.getWindowValue(n, filterLengthSamples, selectedWindow);
                    }
                    double sum = 0.0;
                    for (int i = 0; i < filterLengthSamples; i++) {
                        sum += filterCoefficients[i];
                    }
                    for (int i = 0; i < filterLengthSamples; i++) {
                        filterCoefficients[i] /= sum;
                    }
                }
                try {
                    filteredLeftChannel = maths_used.convolve(leftChannelSamples, filterCoefficients);
                    filteredRightChannel = maths_used.convolve(rightChannelSamples, filterCoefficients);
                    //arr[5] = 10; // This will throw an ArrayIndexOutOfBoundsException

                    System.out.println(("Filtered Left Channel  Data Length :"+filteredLeftChannel.length));
                    System.out.println(("Filtered Right Channel  Data Length :"+filteredRightChannel.length));
                    byte[] newLeftData = maths_used.inttoBytes(filteredLeftChannel,bytesPerSample);
                    byte[] newRightData = maths_used.inttoBytes(filteredRightChannel,bytesPerSample);

                    System.out.println(("Output Left Data Length :"+newLeftData.length));

                    //Interleaving the two Channels
                    newData = maths_used.interleaveData(newLeftData,newRightData);
                    Toast.makeText(MainActivity.this, "Filter Applied", Toast.LENGTH_SHORT).show();
                } catch (Exception e) {
                    Toast.makeText(MainActivity.this, "Try Changing Remez Bands ", Toast.LENGTH_SHORT).show();
                }

            }

        });

        saveButton.setOnClickListener(view -> {
            if (audioData != null) {
                try {
                    // Create the file directory
                    File filepath = new File(getFilesDir(), "MyVolRed");
                    filepath.mkdirs();
// Create the audio file
                    File file = new File(filepath, "audio.wav");
// Get the file URI
                    Uri uri = FileProvider.getUriForFile(MainActivity.this, "com.example.myaudiofilter.fileprovider", file);

                    FileOutputStream fos = new FileOutputStream(file);
                    fos.write(header);
                    fos.write(newData);
                    //Uri uri = FileProvider.getUriForFile(MainActivity.this, "com.example.audio_v2.fileprovider", file);
                    fos.close();// Close the file output stream
                    Toast.makeText(MainActivity.this, "File saved successfully!", Toast.LENGTH_SHORT).show();

                    Intent intent = new Intent(Intent.ACTION_VIEW);
                    intent.setDataAndType(uri, "audio/wav");
                    intent.addFlags(Intent.FLAG_GRANT_READ_URI_PERMISSION);
                    intent.addFlags(Intent.FLAG_ACTIVITY_NEW_TASK);
                    startActivity(intent);

                } catch (IOException e) {
                    Toast.makeText(MainActivity.this, "Failed to save file: " + e.getMessage(), Toast.LENGTH_SHORT).show();
                }
            } else {
                Toast.makeText(MainActivity.this, "No file loaded!", Toast.LENGTH_SHORT).show();
            }
        });
    }

    @SuppressLint("SetTextI18n")
    private final ActivityResultLauncher<Intent> someActivityResultLauncher = registerForActivityResult(
            new ActivityResultContracts.StartActivityForResult(),
            result -> {
                if (result.getResultCode() == RESULT_OK) {
                    Intent data = result.getData();
                    if (data != null) {
                        Uri uri = data.getData();
                        if (uri != null) {
                            try {
                                FileInputStream fileInputStream = (FileInputStream) getContentResolver().openInputStream(uri);
                                audioData = new byte[fileInputStream.available() - 44];
                                fileInputStream.read(header);
                                fileInputStream.read(audioData);
                                numChannels = header[22];
                                bitsPerSample = header[34];
                                bytesPerSample = (int) (header[34] / 8);
                                audioData_Decimal = maths_used.convertByteArrayToIntArray(audioData, bytesPerSample);
                                Toast.makeText(MainActivity.this, "Audio File Loaded Successfully", Toast.LENGTH_SHORT).show();
                                TextView BitsPerSample = findViewById(R.id.BitsperSample);
                                TextView NumChannels = findViewById(R.id.NumChannel);
                                BitsPerSample.setText("Bits per Sample : " + bitsPerSample);
                                BitsPerSample.setTextColor(Color.RED);
                                NumChannels.setText("Number of Channels : " + numChannels);
                                NumChannels.setTextColor(Color.RED);
                                fileInputStream.close();
                            } catch (IOException e) {
                                e.printStackTrace();
                                Toast.makeText(MainActivity.this, "Error: " + e.getMessage(), Toast.LENGTH_SHORT).show();
                            }
                        }
                    }
                }
            });



}
