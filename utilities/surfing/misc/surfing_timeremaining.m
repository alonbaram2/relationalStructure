function [x,y]=surfing_timeremaining(progress,prefix,postfix)
% Prints a status message of the time elapsed and remaining for set of
% operations/
%
% SURFING_TIMEREMAINING(PROGRESS,PREFIX,POSTFIX)
% INPUTS: 
%   PROGRESS    Scalar between 0 and 1 indicating how many of the 
%               operations have been completed (0,1 means none,all)
%   PREFIX   }  Text added before or 
%   POSTFIX  }             after the status message
%
% This function prints the text "took X sec, ETA Y sec", surrounded by
% prefix and postfix if they are given (ETA=Estimated Time of Arrival). If
% prefix is omitted, it is set to the calling function
%
% It is assumed that when the operations were started, the function TIC() 
% has been evoked. The present function calls TOC() to estimate the time 
% remaining using linear extrapolation.
%
% If called with one output argument, S=SURFING_TIMEREMAINING(...), then no
% message is printed but the message is returned in the string S.
%
% If called with two output arguments, [X,Y]=SURFING_TIMEREMAINING(...),
% then the number of seconds it took (X) and the ETA (Y) are returned.
%
% See also TIC, TOC
%
% NNO Dec 2010

if nargin<2, 
    [st,i]=dbstack();
    if numel(st)>=i+1
        prefix=[st(i+1).name ': ']; % caller
    else
        prefix='';
    end
end
if nargin<3, postfix=''; end

tc=toc();
eta=(1-progress)/progress*tc;

msg=sprintf('%stook %d sec, ETA %d sec%s', prefix,round(tc), round(eta),postfix);

switch nargout
    case 0
        fprintf('%s\n', msg);
    case 1
        x=msg;
    case 2
        x=tc;
        y=eta;
end
