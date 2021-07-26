
%% Find the shortest distance between a line and the normal. 
% Reference https://www.mathworks.com/matlabcentral/answers/95608-is-there-a-function-in-matlab-that-calculates-the-shortest-distance-from-a-point-to-a-line

v1 = [0,0,0]; %Must be a 3D vector to work with cross product
v2 = [3,3,0];
pt = [4,1,0];
[xi_bar,x_r,Type, Normal]=point_to_line(pt, v1, v2,"Plot")


function [d, Intercept, Type, m_and_b] = point_to_line(pt, v1, v2,Case)
%%% This function takes in a point (pt) and the endpoints of a line segment
%%% (v1,v2), and returns the shortest distance from the point to line along
%%% the normal as well as the intercept of the normal and line segment. If
%%% the intercept of the normal occcurs outside of the line segment, then 
%%% the function will return the distance to the closest end point. 

%%% NOTE: This function assumes that v1(1)<v2(1)
      tol = 10^-5;
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a); %Find the shortest distance from point to line on the normal
      
      if v1(2)==v2(2)                       %Check if line segment is a horizontal line
          if v1(2)>pt(2)                    %If yes, the intercept can be found easily
            Intercept = [pt(1),v1(2)+d]
          elseif v1(2)<pt(2)
              Intercept = [pt(1),v1(2)-d]
          else
              Intercept = [pt(1),v1(2)]
          end
          
      elseif abs(v1(2)-v2(2))<tol                   %Check if vertical line segment
          "ERROR! Boundary cannot be vertical"
      else
          LineSpecs = polyfit([v1(1),v2(1)],[v1(2),v2(2)],1); %Returns [m,b] for line
          a1=LineSpecs(1);
          b1=LineSpecs(2);
          a2 = -1/a1;
          b2 = (1/a1)*pt(1)+pt(2);
          Intercept = [(-b2+b1)/(a2-a1),a1*((-b2+b1)/(a2-a1))+b1];
      end
      
      if v1(1)<v2(1)                            %Identify which point has the largest and smallest coordinate in x
          MinPnt = v1;
          MaxPnt = v2;
      else
          MinPnt = v1;
          MaxPnt = v2;
      end
      
      if Intercept(1)<MinPnt(1)           %Check if the intercept falls to the left of v1 in the X direction
          Intercept = MinPnt;
          X=[pt;MinPnt];
          d = pdist(X,'euclidean');
      elseif Intercept(1)>MaxPnt(1)       %Check if the intercept falls to the right of v2 in the X direction
          Intercept = MaxPnt;
          X=[pt;MaxPnt];
          d = pdist(X,'euclidean');
      end   
      
      if abs(v1(2)-v2(2))<tol %Check if line segment is horizon)al
          Type = "Horizontal";
          m_and_b = [];
      elseif abs(v1(2)-v2(2))>tol && abs(d)<tol %Check if line segment is angled and pt falls on line segment
          Type = "Angled";
          m_and_b = [a2,b2]
      else                                      %Check if line segment is angled and pt does not fall on line segment
          Type = "Angled";
          m_and_b = polyfit([Intercept(1),pt(1)],[Intercept(2),pt(2)],1); %Returns [m,b] for line
      end
      
      switch Case
          case "Plot"
              plot(pt(1),pt(2), "r*")
              hold on
              plot([v1(1),v2(1)],[v1(2),v2(2)], "Linewidth",2)
              hold on
              plot(Intercept(1),Intercept(2),"b*")
              hold on
              plot([pt(1),Intercept(1)],[pt(2),Intercept(2)], "Linewidth",2)
              xlabel("X")
              ylabel("Y")
              legend ("Point of Interest", "Initial Line Segment", "Intercept between Initial Line Segment and Normal","Line Segment along Normal")
          case "No Plot"
      end
      
      
end



