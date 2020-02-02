classdef traceGen_p2
    methods(Static)
          
        %to reduce the horizontal trace, move the previous point towards
        %different boundary
        %logic indexing
          function [traceXNew,traceYNew] = iterate(interval,traceX,traceY, Vx, Vy, dt)
              T=300;
              [~,numParticle] = size(traceX);
              % assign initial values
              traceXNew = traceX;
              traceYNew = traceY;
              timeArray = linspace(0,dt*interval,interval+1);
              tempArray = zeros(interval);
              LastCollision = zeros(1,numParticle);
              NextCollision = zeros(1,numParticle);
              FreePathHist = zeros(interval,numParticle);
              
              
              %plot init
              figure(2);
              ax1 = subplot(2,1,1);
              ax2 = subplot(2,1,2);
              tag=[];
              %loop over dt
              for i=1:interval
                 
                  [Vx,Vy,LastCollision,NextCollision,FreePathHist(i,:)]=traceGen_p2.scatter(Vx,Vy,T,LastCollision,NextCollision,FreePathHist(i,:));
                  
                  %check whether the partical is going to hit boundary
                  Xnext = traceXNew(i,:)+(Vx*dt);
                  checkX = traceGen_p2.bounceCheck(Xnext,0,200e-9);
                  Ynext = traceYNew(i,:)+(Vy*dt);
                  checkY = traceGen_p2.bounceCheck(Ynext,0,100e-9);
                  
                  %next positions
                  [traceXNew(i,:), traceXNew(i+1,:), Vx ]= traceGen_p2.stepNext(checkX,traceXNew(i,:),Vx, dt,1);
                  [traceYNew(i,:),traceYNew(i+1,:), Vy ]= traceGen_p2.stepNext(checkY,traceYNew(i,:),Vy, dt,0);
                  
                  
%                   plot(traceXNew(i:i+1,:),traceYNew(i:i+1,:));
%                   %remove hold on to see bug in water
%                   hold on;
%                   pause(0.01);
                  
                  
                  color=[1,1,1];
                  for n=1:numParticle
                      %put on ax1 does not work
                    plot(ax2,traceXNew(i:i+1,n),traceYNew(i:i+1,n),'color',color);
                    %remove hold on to see bug in water
                    hold on;
                    pause(0.001);
                    color=color-[0.09,.09,0];
                  end
                  
                  tempArray(i)=traceGen_p2.getTemp(Vx, Vy);
                  plot(ax1, timeArray(1:i),tempArray(1:i));
                  
                  title(ax1,['The average temperature is ',num2str(tempArray(i)),' K'])
                  
                  title(ax2,['The Mean Time of collision is ',num2str(mean(NextCollision-LastCollision)),' s'])
                  delete(tag);
                  tag = annotation('textbox', [0.7, 0.1, 0.1, 0.1], 'String', "MFP: "+num2str(mean(mean(FreePathHist(1:i,:))))+"m");
                  
                  label
              end
          
          end
          
          function checkArray = bounceCheck(pos,limLow,limHigh)
              %mode defines behavior for x,y,etc. 
              tempCheck = zeros(1,numel(pos)) + 1;
              for i=1:numel(pos)
                  element = pos(i);
                  if element<=limLow||element>=limHigh
                      tempCheck(1,i)=-1;
                  end
              end
              checkArray = tempCheck;
              %checkArray is one if particle inside boundary, -1 if across
              %boundary
          end
          
          function [Pos, nextPos, nextVel] = stepNext(checkArray, position, Velocity, dt, mode)
              %mode defines the behavior of the collision
              %THe previous position may be modified for the purpose of
              %jumping over boundaries
              Pos = position;
              switch mode
                  case 0
                      %reflection
                      nextVel = checkArray.*Velocity;
                      nextPos = position + nextVel.*dt;
                      
                  case 1
                      %jump
                      
                      boundary1=0;
                      boundary2=200*1e-9;
                      nextVel = Velocity;
                      nextPos = position + nextVel.*dt;
                      
                      for i=1:numel(nextPos)
                          if nextPos(i)>boundary2 
                              nextPos(i)=nextPos(i)-boundary2;
                              Pos(i)=nextPos(i)-nextVel(i)*dt;
                          end
                          if nextPos(i)<boundary1
                              nextPos(i)=nextPos(i)+boundary2;
                              Pos(i)=nextPos(i)-nextVel(i)*dt;
                          end
                      end
                  otherwise
                      fprintf('matlab NMSL, %d',mode)
              end
              
          end
          function [VxNext, VyNext,LastCollision,NextCollision,FreePath] = scatter(Vx,Vy,Temperature,LastCollision,NextCollision,FreePath)
              %mode defines behavior for x,y,etc. 
              %temperature is for determine the velocity after scatter
              me = 0.26*9.10938215e-31;
              kb = 1.3806504e-23;
              
              dt = 15e-15;
              Tmean = 0.2e-12;
              VThermalMean = sqrt(2*kb*Temperature/me);
              Pscat=1-exp(-dt/Tmean);
              
              for i=1:numel(Vx)
                  if rand()<Pscat
                      AngleParticle = 360*rand(1,1);
                      VThermal = VThermalMean+1e4.*randn(1,1);
                      Vx(i) = VThermal.*cos(AngleParticle);
                      Vy(i) = VThermal.*sin(AngleParticle);
                      LastCollision(i)=NextCollision(i);
                  else
                      NextCollision(i)=NextCollision(i)+dt;
                  end
                  FreePath(i)=(NextCollision(i)-LastCollision(i))*sqrt(Vx(i)^2+Vy(i)^2);
              end
              
              VxNext=Vx;
              VyNext=Vy;
          end
          
          function temp = getTemp(Vx, Vy)
              averageVSq = mean(Vx.^2)+mean(Vy.^2);
              %given average kinetic energy=1/2*mv^2=3/2*kT
              kb = 1.3806504e-23; 
              me = 0.26*9.10938215e-31;
              %use mean of squared velocity, due to the negative v
              temp = 1/kb*(1/2*me*averageVSq);

              
          end
          
    end
end