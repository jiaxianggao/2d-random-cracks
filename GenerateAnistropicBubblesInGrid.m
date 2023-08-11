function [Bubb] = GenerateAnisotropicBubblesInGrid(g,H,min_c,max_c,width_x,width_y,xi,yi,fignr,CONDUCTIVE_BLOB);
%
%
% Generates a bubble distribution
% Aku Seppänen
% 31.1.2000
%
% 

Bubb = ones(length(g),1);
x=g(:,1);
y=g(:,2);

%figure(fignr),PlotSolScaled(g,H,Bubb,0,1),axis image,drawnow;
%button = 1;

%keyboard

for ii=1:length(xi)

  %z = exp(-beta_x(ii)^2*((x-xi(ii)).^2)-beta_y(ii)^2*((y-yi(ii)).^2));
 
  beta_x(ii) = width_x(ii)/6;
  beta_y(ii) = width_y(ii)/6;
z = (1./sqrt(2*pi*beta_x(ii)^2)*exp(-((x-xi(ii)).^2)./(2*beta_x(ii)^2))).*(1./sqrt(2*pi*beta_y(ii)^2)*exp(-((y-yi(ii)).^2)./(2*beta_y(ii)^2)));

  %BuB = 1/max(max(z))*z*rand(1);
  BuB = 1/max(max(z))*z*.5;
  Bubb = Bubb - BuB;
  %figure(fignr),clf,PlotSolScaled(g,H,Bubb,min(Bubb),max(Bubb)),axis image,drawnow;

end;

Bubb = (Bubb-min(min(Bubb)))/(max(max(Bubb))-min(min(Bubb)));
Bubb = min_c + (max_c-min_c)*Bubb;

if CONDUCTIVE_BLOB
  Bubb = min(min(Bubb))+max(max(Bubb)) - Bubb;
end


%close(fignr);
