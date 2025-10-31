function op = buffer_with_truncate(sig,d)
%BUFFER_WITH_TRUNCATE truncates the input signal to an integer
%multiple of the frame size and then buffers it into frames of length d.

% Truncate signal to integer multiple of frame size
sig = sig(1:floor(length(sig)/d)*d); 
% Buffer truncated signal into frames
op=buffer(sig,d);
end