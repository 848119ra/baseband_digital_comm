% student name: Rahele Ahmadian


clc
close all
clear all

%% First define the general system parameters as follows:
SNR = 0:1:15; % Signal-to-noise ratio vector [dB]
%SNR = 15;
T = 1/40e6; % Symbol time interval [s], the symbol rate is 40 MHz
Rs = 40e6; %symbol rate
r = 2; % Oversampling factor (r samples per pulse)
N_symbols_per_pulse = 10; % Duration of TX/RX filter in symbols, Duration of RRC filter in symbols(Nd)
alpha = 0.20; % Roll-off factor (excess bandwidth)
Fs = r/T; % Sampling frequency
Ts = 1/Fs; % Sampling time interval


% GENERATION OF QAM SYMBOLS
% define the number of symbols to be transmitted
N_symbols = 100000; % Number of symbols
% Create a 16-QAM constellation
% Alphabet size
M = 16;
% generally for different alphabet/constellation sizes ():
qam_axis = -sqrt(M)+1:2:sqrt(M)-1;
% generation of a complex constellation:
alphabet = bsxfun(@plus,qam_axis',1j*qam_axis); %help bsxfun
alphabet = alphabet(:).'; % alphabet symbols as a row vector

% Scaling the constellation, so that the mean power of a transmitted symbol
% is one (e.g., with QPSK this is 1/sqrt(2), and for 16-QAM 1/sqrt(10))
alphabet_scaling_factor = 1/sqrt(mean(abs(alphabet).^2));
alphabet = alphabet*alphabet_scaling_factor;

% Number of bits, defined for a fixed alphabet size and number of symbols
N_bits = log2(M)*N_symbols;

% Random bit sequence
bits = randi(2,N_bits,1)-1;





%% Transmitter
cs = 1;
if cs ==1
    nb = 7;

% Generetor Matrix for (7,4) hamming code
P = [ 1 0 1; 1 1 1; 1 1 0; 0 1 1];        
G = [ eye(4) P ];
H = [ P' eye(3) ];

bin = reshape(bits,[],4);

coded = zeros(50000,7);
for j=1:1:length(bits)/4  
    coded(j,:) = mod(bin(j,:)*G,2); 
 
end
coded = reshape(coded,[],1);


elseif cs ==0
    coded=bits;
end

B = reshape(coded,log2(M),[]); %number of bits in each column
q = 2.^(log2(M)-1:-1:0); 
symbol_indices = q*B;

% Gray coded symbol indices.
[Gray_symbol_indices, ~] = bin2gray(symbol_indices, 'qam', M);
% Symbols from the alphabet, based on the Gray-coded symbol indices above.
symbols = alphabet(Gray_symbol_indices+1);

% Plot the symbols
figure(1);
plot(symbols,'ro', 'MarkerFaceColor','r')
axis equal
xlabel('Re')
ylabel('Im')
title('Transmitted symbols')

% Implement the transit filter: Root-Raised-Cosine (RRC) and plot the pulse shape
% Filter generation
gt = rcosdesign(alpha,N_symbols_per_pulse,r,'sqrt');

% Plot the pulse shape of the transmit/receive filter
figure(2);
plot(-N_symbols_per_pulse*r/2*Ts:Ts:N_symbols_per_pulse*r/2*Ts,gt,'b')
hold on
stem(-N_symbols_per_pulse*r/2*Ts:T:N_symbols_per_pulse*r/2*Ts,gt(1:r:end),'ro')
xlabel('time [s]')
ylabel('Amplitude')
title('Transmit/receive RRC filter (pulse shape)')
legend('Pulse shape','Ideal symbol-sampling locations')

%Filter delay correction
symbols_upsampled = zeros(size(1:r*N_symbols*nb/4)); %zero vector initialized for up-sampled sequence
symbols_upsampled(1:r: r*N_symbols*nb/4)= symbols; %symbol insertion

st = filter(gt,1,symbols_upsampled); % Transmitter filtering
st = st(1+(length(gt)-1)/2:end); % Filter delay correction

% Plot the eye-diagram of the generated signal
figure(3); 
hold on 
ks = length(st(1:2000));
for i = 1:2*r:ks % gaussian noise = 2 in 2*r  
     plot(real(st(i:i+2*r)));    
     title('Eye Diagram of the generated signal');    
end 
hold off 
grid on


%figure(3);
%hold on
%for i = 1:2*r:(length(st)-2*r)
% plot(real(st(i:i+2*r)));
%end
%hold off
%grid on

% Plot the transmit signal s(t) in frequency domain
NFFT = 2^14; %FFT size
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts); %frequency vector

% Plot the transmit signal in frequency domain
figure(4);
subplot(211);
plot(f/1e6, fftshift(abs(fft(st, NFFT))));
xlabel('Frequency [MHz]')
ylabel('Amplitude ')
title('TX signal s(t)')
ylim([0 500]);

%% Transmission over channel
% Channel model #6, the channel impulse response is given
h = [ 0.5718 - 0.2667i, -0.4363 - 0.6420i, 0.5645 + 0.1900i, 0.0529 + 0.4435i, -0.1266 - 0.2898i];
% From exercise3, not sure about parameters of FIR filter
h = fir1(4,0.5);
L1 = 1; L2 = 3;
% Channel maximum tap is the first one
%figure
%stem(-L1:L2, abs(h), 'r');
%legend('ISI Channel');
%title('Absolute values of impulse responses');
 

% plot f response of channel 
figure(5);
subplot 211;
% amplitude response of the channel model 6
plot(f/1e6,20*log10(fftshift(abs(fft(h,NFFT)))));    
xlabel('Frequency [MHz]'); 
ylabel('Amplitude Response [dB]');  
legend('Channel Model 6'); 

figure(5);
subplot 212;
%phase response of the channel model 6
plot(f/1e6,fftshift((angle(fft(h,NFFT)))));   
xlabel('Frequency [Hz]');  
ylabel('Phase Response [degrees]'); 
legend('Channe Model 6');

% Pass signal through channel

st_fading = conv(st, h); % convolution of transmitted signal and channel impulse response

% Generate the noise vector
% Complex white Gaussian random noise
n = (1/sqrt(2))*(randn(size(st_fading)) + 1j*randn(size(st_fading)));
P_s = var(st_fading); % Signal power
P_n = var(n); % Noise power  
% Defining noise scaling factor based on the desired SNR:
noise_scaling_factor = sqrt(P_s/P_n./10.^(SNR./10)*(r/(1+alpha)));

% Initialization for RX signal matrix, where each row represents the
% received signal with a specific SNR value
rt = zeros(length(SNR), length(st_fading));
% Received signal with different SNR values
for ii = 1:1:length(SNR)
 rt(ii,:) = st_fading + noise_scaling_factor(ii)*n;
end

%Received signal with the noise when the SNR is equal to the last value of
%the SNR vector
figure(4);
subplot(212);
plot(f/1e6, fftshift(abs(fft(rt(end,:), NFFT))));
xlabel('Frequency [MHz]')
ylabel('Amplitude')
title(['RX signal r(t) (SNR = ', num2str(SNR(end)), ' dB)'])
ylim([0 500]);

%% Receiver
% Signal filtering, Filter the received signal r(t) with the receive filter (RRC similar to TX) 
ft = gt;

% Initialization for the received symbol matrix, where each row represents
% the symbols with a specific SNR value
qk = zeros(length(SNR), N_symbols - N_symbols_per_pulse);

% Filtering and sampling
for ii = 1:1:length(SNR)
 qt = filter(ft,1,rt(ii,:)); % Receiver filtering 
 qt = qt(1+(length(ft)-1)/2:end); % Filter delay correction
 
end

% Plot the eye-diagram of the received signal
figure(6);
hold on; 
for i=1:50*r:(length(qt)-2*r) 
plot(real(qt(i:i+2*r)),'b');  
end
hold off;
title('Eye Diagram RX and Noise signals)');   
xlabel('Relative Time [s]');    
ylabel('Amplitude'); 

% Downsampling to the symbole rate, remove oversampling
rt_downsampled=qt(1:r:end);

figure(7);
hold on; 
plot(rt_downsampled,'ob'); 
xlabel('Real');
ylabel('Imagianry');   
legend('Rx Constellation','Constellation 16_ QAM ');
hold off


%Channel Estimation:

M_ch = 30; %number of reference symbols used for channel estimation; 
estimate_length = 7; %how long is the channel estimate's impulse response
A_conv=convmtx(st(2:M_ch+1).',estimate_length); % Convolution matrix
p_LS=((A_conv'*A_conv)\(A_conv'))*rt(1:size(A_conv,1)).'; % LS solution

figure(8); 
hold on; 
[H,f]=freqz(p_LS,1,-Rs/2:Rs/200:Rs/2,Rs); 
plot(f/1e6,20*log10(abs(H)),'r'); 
legend('Channel','LS Channel Estimate');
title('Channel Estimation');

figure(9);
hold on; 
stem(-L1:length(h)-L1-1,abs(h),'k'); 
stem(-L1:length(p_LS)-L1-1,abs(p_LS),'r'); 
legend('Channel','LS channel estimate'); 
title('Absolute values of the impulse responses')


% Channel equalizer using LMS algorithm
beta = 0.001; % step-size of the algorithm
c_LMS = zeros(31,1); % equalizer coefficients, initializations
for i = 16:length(rt_downsampled)-15 
 rk = flipud(rt_downsampled(i-15:i+15).'); % Received signal vector
 Ek(i) = st(i) - c_LMS.'*rk; % Error signal, we assume a known symbol sequence
 c_LMS = c_LMS + beta*Ek(i)*conj(rk); % LMS update !
end

Rx_eq = filter(c_LMS,1,rt_downsampled);

Rx_eq = Rx_eq(1+(length(c_LMS)-1)/2:end);

figure(10);
hold on; 
plot(abs(Ek)); 
title('Convergence behavior of the LMS-algorithm.'); 
ylabel('LMS error'); 
xlabel('Iteration index'); 

figure(11); 
hold on; 
stem(abs(conv(h,c_LMS))); 
title('Effective impulse response (abs) of the equalized system ') 

figure(12); 
hold on; 
[H,f]=freqz(c_LMS,1,-Rs/2:Rs/200:Rs/2,Rs); 
plot(f/1e6,20*log10(abs(H)),'r'); 
[H,f]=freqz(conv(c_LMS,h),1,-Rs/2:Rs/200:Rs/2,Fs); 
plot(f/1e6,20*log10(abs(H)),'g'); 
legend('Channel','LMS Equalizer','Total Response (LMS)'); 
title('Channel vs LMS Equalizer');

figure(13); 
hold on; 
plot(filter(c_LMS,1,rt_downsampled),'rx'); % Equalized constellation
%plot(symbols,'.k', 'MarkerFaceColor','r') % Originally generated symbols constellation
axis equal
xlabel('Re')
ylabel('Im')
title('Transmitted symbols')
plot(rt_downsampled,'ob'); % Received signal constellation after channel 
xlabel('Real');
ylabel('Imagianry');   
legend('LMS Equalized constellation','Constellation 16_ QAM ','Rx Constellation');
hold off

% BER calculation
% Initialization


%% 2.4 Error control coding
cs = 1;
if cs ==1
    nb = 7;

P = [ 1 0 1; 1 1 1; 1 1 0; 0 1 1];        
G = [ eye(4) P ];
H = [ P' eye(3) ];

bin = reshape(bits,[],4);

coded = zeros(50000,7);
for j=1:1:length(bits)/4  
    coded(j,:) = mod(bin(j,:)*G,2); 
 
end
coded = reshape(coded,[],1);


elseif cs ==0
    coded=bits;
end

B = reshape(coded,log2(M),[]); %number of bits in each column
q = 2.^(log2(M)-1:-1:0); 
symbol_indices = q*B;

[Gray_symbol_indices, ~] = bin2gray(symbol_indices, 'qam', M);
symbols = alphabet(Gray_symbol_indices+1);



alphabet_error_matrix = abs(bsxfun(@minus,alphabet.',Rx_eq));
[~,estimated_Gray_symbol_ind] = min(alphabet_error_matrix);
estimated_symbol_indices = gray2bin(estimated_Gray_symbol_ind-1,'QAM',16);
estimated_bit_blocks = rem(floor((estimated_symbol_indices(:))*2.^(1-log2(M):0)),2)';
estimated_bits = estimated_bit_blocks(:);

comp = zeros(length(coded),1);
dif = zeros((length(coded)-length(estimated_bits)),1);
comp(:,1)= [ estimated_bits(:); dif ];
ecic = reshape(comp,[],7);
estimated_bits = comp;

bit_errors = estimated_bits ~= coded(1:length(estimated_bits)); 
         
   
   BER = mean(bit_errors); 
n=7;

if cs == 1
syb = zeros(((length(coded))/n),3);
for j=1:1:length(coded)/n
    
    syb(j,:) = mod(ecic(j,:) * H',2);
end

find = 0;
errc = 0;  
epci = zeros(1,length(estimated_bits));
epri = zeros(1,length(estimated_bits));
end

for i = 1:1:length(estimated_bits)/n
    
    for j=1:1:n
       erv= zeros(1,n);
        erv(j) = 1;
        sch = mod( erv * H',2);
        
        
        if sch == syb(i,:)
           
            
            errc = errc+1; 
            
            epci(errc) = i; 
            epri(errc) =j;
        end
       
        end
    end  
 
  
    ccode = ecic ;

    for j=1:1:errc
    
        ccode( epci(j),epri(j) ) = not( ecic( epci(j),epri(j) ));
    
    end
 
    decoded_bits = reshape(ccode(:,1:4),[],1);

    if cs == 0
     
        ccode = estimated_bits ;
        decoded_bits =reshape(ccode,[],1);    
    end 
    %% 
bits = decoded_bits ~= coded(1:length(decoded_bits));
    BER = mean(bits);



% Compare the simulated and theoretical results.
figure(14);
semilogy(SNR, BER, 'LineWidth', 2);
hold on;
title('Bit error rate')
xlabel('SNR [dB]')
ylabel('BER')
legend('Simulated BER with Gray coding',...
 'Theoretical bit error probability (approx.)','Location', 'SouthWest');


% Question1: What is the effect of the roll-off factor on the transmitted signal in general?
% By increasing the roll-off factor, from 0.2 to 1, the excess bandwidth of the transmitted signal will be increased.
% For roll-off factor = 0.2, the excess bandwidth is 20% and for roll-off factor = 1 it will be increased to 100%.


% Question2: Why the eye-diagram of the transmitted signal doesn%t show pure/perfect eye opening?
% In my simulation parameters, the roll-off factor is set to 0.2, and when I change it to 1, then the eye-diagram 
% had perfect eye opening. So by roll-off factor of 0.2, 20% excess bandwidth, the pulse shape is more similar to 
% sinc shape and more sensitive to timing errors, so the horizontal openings of eye-diagram are smaller and not perfect.

% Question3: Is it feasible to use a fixed impulse response as channel model for a radio channel in 
% mobile communication systems? Why? Why not?
% No, it is not feasible to use a fixed impulse response as channel model for mobile communications,
% because of multipath delays. As the channel is dynamic and there will be multipath delays and the impulse response 
% of the channel changes and using a fixed impulse response is impractical. Here by using FIR filter for the channel, 
% the effects of multipath delays on the received signal has been mitigated.

% Question4: What is the meaning of the eye-diagram? What can be observed in an eye-diagram?
% An eye-diagram consists of several overlapping symbol level synchronized waveform traces. It is assumed that symbols
% are random and independent, in which case all possible symbol combinations are present. ISI and sensitivity to timing
% synchronization errors can be seen in the eye-diagram. If the vertical eye opening in the eye-diagram is greater, 
% The system is more immune to noise. If the horizontal eye opening in the eye-diagram is smaller then the system is more
% sensitive to timing synchronization. By increasing the roll-off factor, the horizontal eye opening of the eye-diagram
% will be increased. In my simulation, when I plot the eye diagram of the received signal, there is no horizontal and
% vertical eye opening. And it shows that the system is sensitive for timing synchronization and, more vulnerable to noise.

% Question5: What is the idea of matched filtering in general? Consider true matched filter
% (matched to the received pulse) versus matched filter to the transmitted pulse?
% Matched filtering is a technique used in digital communication systems for efficient detection and recovery of the 
% transmitted signals in the presence of noise. The matched filter is used in both transmitter and receiver side. 
% To shape and reshape the signal for better SNR in the receiver, and improving the overall system performance. 
% A true matched filter refers to a filter that perfectly matches the shape and characteristics of the transmitted
% pulse. Implementing true matched filter can be challenging in practice because of channel impairments and limitation
% in hardware. The matched filter to the transmitted pulse is a filter that approximates the true matched filter based 
% on the knowledge of the transmitted pulse shape, and it is designed in a way that maximizes the SNR of the receiver. 

% Question6: What is the idea of an equalizer in a digital communication system?
% Channel equalizer is an alternative to handle ISI in the receivers. The idea is to reduce or remove the amount
% of ISI in the received discrete time signal using a discrete time linear filter (equalizer), after this make decision
% about the transmitted symbols. Practical equalizer filters have finite complexity, so it may not always be feasible 
% to completely remove all ISI. The equalizer%s primary function is to estimate and reverse the effects of the channel
% on the received signal.

% Question7: What performance would you expect over the given SNR range using a soft-decision decoder instead 
% of hard-decision decoder. Why, what is the trade-off?
% If we use soft-decision decoder instead of hard-decision decoder (based on Hamming code H(7,4)) that was used here,
% the BER will be decreased and error control coding performance will be better, as soft decoders are better in
% performance than hard decoders, but as in soft decoding symbol detection and channel decoding are combined, 
% soft decoding is more complicated than hard decoding and that is the trade-off here, between performance and
% complexity.














