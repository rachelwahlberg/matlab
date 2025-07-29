function playLocalAudioFile(bellSoundPath)
    % Check if the file exists at the specified path
    if exist(bellSoundPath, 'file') == 2
        % Read the audio data from the file
        [y, Fs] = audioread(bellSoundPath);
        
        % Play the audio
        sound(y, Fs);
    else
        % If the file does not exist, display an error message
        fprintf('Error: File does not exist at the specified path.\n');
    end
end