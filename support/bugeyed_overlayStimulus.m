function im = bugeyed_overlayStimulus(im, visField, stimPara, t)
   
    switch stimPara.type
        case 'MovingDot'
            diam = stimPara.initialDiameter;
            pos = stimPara.initialPosition + t*stimPara.movementSpeed;
    
        case 'LoomingDot'
            pos = stimPara.initialPosition;
            d   = max([0 1 - t/stimPara.timeToContact]); % distance of the virtual object (0-1)
            diam = 2 * atand(tand(stimPara.initialDiameter/2)/d);

        otherwise
            error('Unknown stimulus type: %s', stimPara.type)
    end
    im = bugeyed_overlayDot(im, visField, pos, diam, stimPara.opacity);

end