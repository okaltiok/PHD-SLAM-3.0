% This code is based on an implementation posted by Edric Ellis on 10 Jun 2019 at:
% "https://se.mathworks.com/matlabcentral/answers/465911-parfor-waitbar-how-to-do-this-more-cleanly"
% and slightly modified for the desired purpose.

classdef PoolWaitbar < handle
    properties (SetAccess = immutable, GetAccess = private)
        Queue
        N
        start_time
    end
    properties (Access = private, Transient)
        ClientHandle = []
        Count = 0
    end
    properties (SetAccess = immutable, GetAccess = private, Transient)
        Listener = []
    end
    methods (Access = private)
        function localIncrement(obj)
            t = now;
            obj.Count = 1 + obj.Count;
            waitbar(obj.Count / obj.N, obj.ClientHandle);
            
            str = [...
                'Round: ',num2str(obj.Count),' / ', num2str(obj.N), ', Time: ', ...
                sprintf('%02d',round(hour(t-obj.start_time))),':', ...
                sprintf('%02d',round(minute(t-obj.start_time))),':', ...
                sprintf('%02d',round(second(t-obj.start_time))) ...
                ];
            
            set(obj.ClientHandle.Children(1).Title,'String',str)
            pause(0.1)
        end
    end
    methods
        function obj = PoolWaitbar(N, message)
            if nargin < 2
                message = 'Simulation in progress, please wait ...';
            end
            obj.N = N;
            obj.start_time = now;
            obj.ClientHandle = waitbar(0, message);
            obj.Queue = parallel.pool.DataQueue;
            obj.Listener = afterEach(obj.Queue, @(~) localIncrement(obj));
        end
        function increment(obj)
            send(obj.Queue, true);
        end
        function delete(obj)
            delete(obj.ClientHandle);
            delete(obj.Queue);
        end
    end
end