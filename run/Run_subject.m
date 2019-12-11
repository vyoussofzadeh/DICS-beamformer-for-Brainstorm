switch task
    case 1
        % - Auditory definition naming
        tag = 'DFN'; tag1 = 'dfn';
        Evnt_IDs = 1; % questions
        disp('1: 2019')
        disp('2: 2018');
        disp('3: 2017');
        disp('4: 2016')
        disp('5: 2015');
        disp('6: 2014');
        disp('7: 2013')
        disp('8: 2012');
        disp('9: 2011');
        disp('10: all');
        year = input('Year data were acquired: ');
    case 2
        % - Visual picture naming
        tag = 'PN';
        Evnt_IDs = 3; % 3: images, 2: scrambled images
        disp('1: 2019')
        disp('2: 2018');
        disp('3: 2017');
        year = input('Year data were acquired: ');
end

%%
clear ytag;
switch year
    case 1
        stag = {'19'};ytag = {'19'};
    case 2
        stag = {'18'};ytag = {'18'};
    case 3
        stag = {'17'};ytag = {'17'};
    case 4
        stag = {'16'};ytag = {'16'};
    case {5,6}
        stag = {'up'};ytag = {'14_5'};
        %     case 6
        %         stag = {'up'};ytag = {'14'};
    case 7
        stag = {'13'};ytag = {'13'};
    case 8
        stag = {'12'};ytag = {'12'};
    case 9
        stag = {'11'};ytag = {'11'};
    case 10
        stag = {'all'};ytag = {'all'};
end
disp('============');