;function y = string_fdtd_s1889125(opts,phys_param,sim_param)

%-------------------------------------------------------------------------%
%_*Program Description*_ 
%Designing a function that carries out an FDTD simulation of a stiff string, 
%which takes in a set of options and parameters and returns a monophonic 
%output audio signal.This function can be used to generate a single string
%sound with the help of struct for parameter input.

%-------------------------------------------------------------------------%
   %print options and parameters
   opts;phys_param;sim_param;

   %copy over parameters, taking into account some options
   T = phys_param.T;          % tension (N)
   r = phys_param.r;          % string radius (m)
   if opts.add_stiffness
      E = phys_param.E;       % Young's modulus (Pa)
   else
      E = 0;
   end
   rho = phys_param.rho;      % density (ρ) (kg/m^3) 

   T60 = phys_param.T60;      % T60 (s)
   L = phys_param.L;          % length (m)

   SR = sim_param.SR;         % sample rate (Hz)

   Tf = sim_param.Tf;         % duration of simulation (s)

   xi = sim_param.xi;         % coordinate of excitation (normalised, 0-1)
   famp = sim_param.famp;     % peak amplitude of excitation (N)
   dur = sim_param.dur;       % duration of excitation (s)
   exc_st = sim_param.exc_st; % start time of excitation (s)
   xo = sim_param.xo;         % coordinate of output (normalised, 0-1)
%-------------------------------------------------------------------------%
  %                         Derived parameters 
%-------------------------------------------------------------------------% 
   
% string cross-sectional area
   A = pi*r^2;  
% string moment of intertia
  I = 0.25*pi*r^4;                  
% wave speed
  c = sqrt(T/(rho*A));
% loss parameter (σ)     
   sig = 6*log(10)/T60; 
% stiffness constant (κ)
   kappa = sqrt(E*I/(rho*A)); 
% time step         
   k = 1/SR;                        
% minimal grid spacing for stability
   hmin = sqrt(0.5* (c.^2*k^2+sqrt(c.^4*k^4+16*kappa.^2.*k.^2)) );
% number of grid points to update
   N = floor(L/hmin) ;
% actual grid spacing used
   h = L/(N);             
  
   assert(h>=hmin)                     %for stability
   assert(sig>=0)                      %for stability

 % Courant number (λ)
   lambda = c*k/h;   
 % numerical stiffness constant (μ)
   mu = (k*kappa/h^2); 
 %Grid points update
   N=N-1;
   
%-------------------------------------------------------------------------%
%                                I/O 
%-------------------------------------------------------------------------%
% number of time steps for simulation (based on Tf, SR)   
    Nf = floor(Tf*SR);  
% grid index of excitation (based on xi,N,L)
   li =floor(xi.*N/L)+1;  
% grid index of output (based on xo,N,L)
   lo = floor(xo.*N/L)+1; 
   
%-------------------------------------------------------------------------%
%                       Error check and update
%-------------------------------------------------------------------------%
   
if any([(xi*L)<h (L-(xi*L))<h])
  warning([sprintf('xo is too close to string boundary\n')]);
  
  xi = abs(round(xi)-h);
  if li<=2
    li = 3;
  elseif li>=N
    li = N-1;
  end
end

   assert(is_pinteger(li))
   assert(is_pinteger(lo))
%-------------------------------------------------------------------------%
%                        create force signal
%-------------------------------------------------------------------------%
   
% input force signal
   f = zeros(Nf,1);
% duration of force signal, in samples
    durint = floor(dur*SR);     
% start time index for excitation
    exc_st_int = (floor(exc_st*SR))+1; 
% sample values of force
    durf = exc_st_int:exc_st_int+durint-1;  
 % force foefficient    
    f_coeff = (k^2./(h*rho*A));         
% force vector
% FOR 'STRUCK' full hann window
if strcmp(opts.input_type,'struck') 
for zz=exc_st_int:exc_st_int+durint-1
f(zz) = famp.*0.5.*(1-cos((2)*pi.*(zz./durint))).*f_coeff;
end
% FOR 'PLUCKED', half hann window
   elseif strcmp(opts.input_type,'plucked')
       for zz=exc_st_int:exc_st_int+durint-1
   f(zz) = famp.*0.5.*(1-cos(pi.*(zz./durint))).*f_coeff;
       end
end

%-------------------------------------------------------------------------%
%                        time/state variables
%-------------------------------------------------------------------------%
% state at time index n+1   
    u0 = zeros(N,1); 
% state at time index n
   u1 = zeros(N,1); 
% state at time index n-1
   u2 = zeros(N,1);        
% output vector
   y = zeros(Nf,1);           
%start and end l-index for for-loop update
   lstart = 3; 
   lend = N-2; 

   assert(is_pinteger(lstart)) 
   assert(is_pinteger(lend)) 
%-------------------------------------------------------------------------%
   %                     main loop and calculation
%-------------------------------------------------------------------------%
   tic;
 %Coefficients pre computation
    cf1= 1+sig*k;
    cf2 = 2-2*lambda^2-(mu^2)*6 / cf1;
    cf3 = 1-sig*k /cf1;
    cf4 = (lambda^2)+(mu^2)*4 /cf1;
    cf5 = mu^2 /cf1;
   for n=1:Nf
  %Using for loop
      % interior update
      if opts.useforloop
         for l = lstart:lend
                u0(l) = (u1(l)*cf2 - u2(l)*cf3 + cf4*(u1(l+1) + u1(l-1))-...
                        cf5*(u1(l+2) + u1(l-2)) );
          end    
  %vectorized
      else 
           u0(lstart:lend) = (u1(lstart:lend)*cf2 - u2(lstart:lend)*cf3 + cf4*(u1(lstart+1:lend+1) + u1(lstart-1:lend-1))-...
                cf5*(u1(lstart+2:lend+2) + u1(lstart-2:lend-2)) );
      end
      
      
      
%-------------------------------------------------------------------------%
%                           boundary updates
%-------------------------------------------------------------------------%
      if strcmp(opts.bctype,'clamped')
          u0(1)=cf1\(2*u1(1)-u2(1)*(1-sig*k) +lambda^2*(-2*u1(1)+u1(2))...
              -mu^2*(u1(3)-4*u1(2)+6*u1(1)));
          u0(2) = cf1\(2*u1(2)-u2(2)*(1-sig*k)+lambda^2*(u1(1)-2*u1(2)+u1(3))...
              -mu^2*(u1(4)-4*u1(3)+6*u1(2)-4*u1(1)));
          

            u0(N)=0;
       
      elseif strcmp(opts.bctype,'simply_supported')
          u0(1) = cf1\ (2*u1(1)-u2(1)*(1-sig*k)+lambda^2*(-2*u1(1)+u1(2))...
              -mu^2*(u1(3)-4*u1(2)+5*u1(1)));
          u0(2) = cf1\ (2*u1(2)+lambda^2*(u1(1)-2*u1(2)+u1(3))...
              -mu^2*(u1(4)-4*u1(3)+6*u1(2)-4*u1(1))-u2(2)*(1-sig*k));
          u0(N-1) = cf1\ (2*u1(N-1)+lambda^2*(u1(N-2)-2*u1(N-1)+u1(N))...
              -mu^2*(u1(N-3)-4*u1(N-2)+6*u1(N-1)-4*u1(N))-u2(N-1)*(1-sig*k));
          u0(N) = cf1\ (2*u1(N)-u2(N)*(1-sig*k)+lambda^2*(u1(N-1)-u1(N)*2)-...
              mu^2*(-5*u1(N-1)+6*u1(N)+u1(N-2)));
         
      end
  % send in input
      u0(li) = u0(li)+f(n);

  % read output
      if strcmp(opts.output_type,'displacement')
         y(n) = u0(lo);
      elseif strcmp(opts.output_type,'velocity')
          
 % 'velocity'' read-out 
         y(n) = (u0(lo)-u1(lo))*SR; 
      end
%-------------------------------------------------------------------------%
%                               plotting
%-------------------------------------------------------------------------%

      if (opts.plot_on)
         % draw current state
         if n==1
            figure
            h1=plot([1:N].'*h, u0, 'k');
            axis([0 L -0.005 0.005])
            xlabel('position (m)')
         else
            set(h1,'ydata',u0)
            drawnow
         end
         fprintf('n=%d out of %d\n',n,Nf);
      end

      % shift states to step forward in time
      u2 = u1;
      u1 = u0;
   end

   %read last samples of output
   for n=Nf-4:Nf
      fprintf('y(%d) = %.15g\n',n,y(n));
   end
   toc

   %%%%% plot spectrum
   if (opts.plot_on)
      figure
      yfft = 10*log10(abs(fft(y)));
      plot([0:Nf-1]'/Nf*SR, yfft, 'k')
      xlabel('freq (Hz)')
   end
end

%is positive integer?
function y=is_pinteger(x)
   y=((mod(x,1)==0) && (x>0));
end
