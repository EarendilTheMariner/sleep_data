function R = SimulateNetwork(W,t_steps,varargin)
    f_V_func = @tanh_f_V;
    tau_V = ones(size(W,1),t_steps)*50;% Slower Excitation
    %tau_V(sum(W,1) < 0,:) = 20; % Faster Inhibition
    noise_ampl = 0.4;
    seed = 5;
    I_e = ones(size(W,1),t_steps)*20;
    V_init = 0;
    threshold = 20;
    gain = ones(size(W,1),t_steps).*1;
    gainfeed = 0.2;
    f_max = 50;
    dfeed = 15;
    Source = ones(size(W,1));
    Target = zeros(size(W,1));
    Position = ones(size(W,1));
    Save = 0;
    dt = 10;

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'tau_V'
                tau_V = varargin{ii+1};
            case 'noise_ampl'
                noise_ampl = varargin{ii+1};
            case 'seed'
                seed = varargin{ii+1};
            case 'I_e'
                I_e = varargin{ii+1};
            case 'V_init'
                V_init = varargin{ii+1};
            case 'threshold'
                threshold = varargin{ii+1};
            case 'gain'
                gain = varargin{ii+1};
            case 'fmax'
                f_max = varargin{ii+1};
            case 'Position'
                Position = varargin{ii+1};
            case 'ClassSource'
                Source = varargin{ii+1};
                LeftSource = Source & Position;
                RightSource = Source & ~Position;
            case 'ClassTarget'
                Target = varargin{ii+1};
                LeftTarget = Target & Position;
                RightTarget = Target & ~Position;
            case 'gainfeed'
                gainfeed = varargin{ii+1};    
        end
    end

    rng(seed);
    N = size(W,1);
    R = zeros(t_steps, N);  
    V = zeros(N, t_steps);
    V(:, 1) = V_init;
    R(1, :) = V_init;
    I_noise = normrnd(0.,noise_ampl,[1 (t_steps)*N]);

    tic
    for t = 2:t_steps
       for i = 1:N
            if t > dfeed && Target(i)
                if(LeftTarget(i))
                    feed = gainfeed(i)*mean(R(t-dfeed,LeftSource));
                elseif(RightTarget(i))
                    feed = gainfeed(i)*mean(R(t-dfeed,RightSource));
                end
            else
                feed = 0;
            end
            I_rec = dot(W(i,:), R(t-1,:));
            I_tot = I_rec + I_noise(t*i) + I_e(i, t) + feed;
            dV = (1. / tau_V(i, t)) * (-V(i,t-1) + I_tot);
            V(i,t) = V(i,t-1) + dV;
            R(t, i) = f_V_func(V(i,t), threshold, gain(i, t), f_max);
       end      
    end
    toc
end

%% Helper Functions 
function val = tanh_f_V(V,threshold,gain,fmax) %20,1,100
    if (V-threshold)<=0.
        f = threshold*tanh(gain*(V-threshold)./threshold);
    elseif (V-threshold)>0
        f = fmax*tanh(gain*(V-threshold)./fmax);      
    end
    val = f+threshold;
end