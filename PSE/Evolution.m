% The form of the FTCS 

function pddt = Evolution(pd)
    global U V alpha beta numPoints_phi_1 numPoints_theta_1 dphi_1
    pd_U = U.*pd;
    pd_V = V.*pd;
    pd_copy = pd;
    for x = [1:numPoints_phi_1+1]
        N = round(x+(pi)/dphi_1);
        if N>(numPoints_phi_1+1)
            N = N-(numPoints_phi_1+1);
        end

%%% Find a way to solve the integral over a non-analytic point. Maybe interpolation, or something else
%         if  abs(round((pd_copy(1,x)-pd_copy(1,N)),3))>0
%             pd_copy(1,x) = 1/2*(pd_copy(1,x)+pd_copy(1,N));
%         end
%         if abs(round((pd_copy(numPoints_theta+1,x)+pd_copy(numPoints_theta+1,N)),2))>0
%             pd_copy(numPoints_theta+1,x) = 1/2*(pd_copy(numPoints_theta+1,x)+pd_copy(numPoints_theta+1,N));
%         end


        pd(1,x) = pd_copy(1,x)+beta/2*(pd_V(2,x)-pd_V(2,N));
        pd(numPoints_theta_1+1,x) = pd_copy(numPoints_theta_1+1,x)+beta/2*(pd_V(numPoints_theta_1,N)-pd_V(numPoints_theta_1,x));
    end
    for x = [2:numPoints_theta_1]

%%%  Find a way to solve the integral over a non-analytic point. Maybe interpolation, or something else
%         if abs(round((pd_copy(x,1)-pd_copy(x,numPoints_phi+1)),2))>0
%             pd_copy(x,1) = 1/2*(pd_copy(x,1)+pd_copy(x,numPoints_phi+1));
%             pd_copy(x,numPoints_phi+1) = 1/2*(pd_copy(x,2)+pd_copy(x,numPoints_phi));
%         end

        pd(x,1) = pd_copy(x,1) + alpha/2*(pd_U(x,2)-pd_U(x,numPoints_phi_1))+beta/2*(pd_V(x-1,1)-pd_V(x+1,numPoints_phi_1));
        pd(x,numPoints_phi_1+1) = pd(x,1);
        
    end
    for x = [2:numPoints_theta_1]
        for y = [2:numPoints_phi_1]   
            pd(x,y) = pd_copy(x,y) + alpha*(pd_U(x,y+1)-pd_U(x,y-1))/2+beta*(pd_V(x+1,y)-pd_V(x-1,y))/2;
        end
    end

    pddt = pd;
end