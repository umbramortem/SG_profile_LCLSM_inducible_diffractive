
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   El presnete codigo es una parte importante en la optencion de 
%
%
%    resultados del manuscrito titulado "Simulated LCSLM with inductible 
%    Diffractive theory to display super-Gaussian array applying the
%    irradiance transport equation", sometido al Journal Optik/Optics at 
%    Elsevier
%
%    Correspondings Authors:
%    dr.j.a.arriaga.hernandez@gmail.com ;  jesus.arriagahdz@correo.buap.mx 
%    b.cuevas.otahola@gmail.com   ;   bolivia.cuevasotahola@viep.com.mx 
%
%    El script se dividio en 6 secciones como se driscribio en el manuscrito
%    ademas del orden de obtencion de resultados.
%    En la Seccion I obtenemos la base de la LCD (LCSLM).
%    En la Seccion II obtnemos la LCD y la comprobamos obcervando los test
%    chart obtenidos de los siguientes enlaces:
%    https://www.hwalworks.com/The-Year-of-Resolution-Test-Chart
%    https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=4338
%    En la Seccion III obtenemos la LCD considerando el borde de su
%    pantalla, donde se peude simular cualquier tamaño en LCD y pantalla.
%    En la Seccion IV se simulas LCD de 1x1, 2x2, 3x3 pixeles para obtener
%    su patron de difraccion por Fourier Transform.
%    En la Seccion V obtenemos los difrraction pattern para el caso de LCD 
%    de 1x1, 2x2, 3x3 y 512x512 pixeles mediante los calculos del pdf 
%    "Supplementary Document"
%    En la Seccion VI obtenemos la SGR (Super Gaussian Riling) como una 
%    rejilla de franjas con perfil Super-Gaussiano.
%
%
%    Debe tenerse en consideracion que los resultados se obtuvieron en una
%    pc Gamming con procesador i7 de 16 nucleos con 64 gb en RAM y una
%    taregta grafica de con 1280 cuda cores y 4gb GDDR6.
%    Los resultados se obtuvieron para matrices de 1024x1024 por lo que se
%    debe tener en consideracion el poder de computo para emular los
%    resultados. Es facil reescalar los resultados para matrices de 512x512
%    pero se deben considerar los parametros a modificar.
%
%
%
%
%
%
%
%%
%%%% corregir todo por secciones
clear all; close all;



%%   Section I
%%%% simple Difraction Aperture Basis

% Pixel size in millimeters 
SquareSize = 200;
a = SquareSize;
% Matrix size in millimeters 
MatrixSize = 3*SquareSize;
b = MatrixSize;
Window = zeros(b,b);

% square aperture
a1 = (b/2)-(a/2);       a2 = (b/2)+(a/2);
% Difracctive aperture basis
Window(a1:a2,a1:a2) = 1;

figure;
imshow(Window); title('Aperture Diffraction Basis');


%% Section II
%%%%
% Relation to determine the pixeal number of the LCD between Acpix and PT,
% N=AcPix/PT
% Mask plane size in millimeters
TT = 560;       %For the Fig. 5, AcPix =560.
% Only Diffraction mask (2D space with square apertures) of active pixel
AcPix = 512;    %For the Fig. 5, AcPix =512.
T = AcPix;
% Pixel size
PT = 1;
% Square pixel number
N = int16(fix(T/PT));
% Pixel real total number, transparets ans opaque
nn = 1 + (2*N);
% Amplification or magnification factor
nnn = 4;
% Visualmente, que tanta resolucion se necesita, es decir, el tamaño de
% pixel real contra un tamaño que visualmente sea aceptable en un grafico.
% La interpretacion de la mascara TxT en una matriz de tamaño MxM. Donde T
% es el tamaño del pixel considerando que matlab en la matrix coloca el
% numero 1 a cualquier valor de entrada, ie, T = 1*n, con n entero tamaño
% visual del pixel en el grafico
T = nnn*nn;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Bouilding just de difraction mask
msk = zeros(T,T);
% first step, row basis
msk(1:nnn,:) = 0;
for i=1:N
    ii = ((2*i) -1)*nnn;
    jj = ii + nnn;
    msk( (nnn+1):(2*nnn),(ii + 1):jj ) = 1;
end

maskbasis(1:nnn,1:T) = 0;
maskbasis( (nnn+1):(2*nnn),1:T) = msk( (nnn+1):(2*nnn),:);
for i=1:N
    ii = (i - 1)*2*nnn + 1;
    
    msk( ii:(2*i*nnn), 1:T ) = maskbasis;
end
% Only Diffractive Mask 
DiffMask = msk;
figure;
imshow(msk); title('Diffraction aquare apertures mask Only');




% Simple example of displaying an image in shades of gray tones trought the
% LCD simulated
% Load image file, test charts
imag1 = imread('C:\Users\jesus\OneDrive\Escritorio\Trabajo ETI Super gaussian\Figures\LCD Image Test\ImageLCDTest1.png');

imag2 = imread('C:\Users\jesus\OneDrive\Escritorio\Trabajo ETI Super gaussian\Figures\LCD Image Test\ImageLCDTest2.png');

% Image resize to square matrix TxT
Imag1 = (imresize(imag1, [T T]));
Imag2 = (imresize(imag2, [T T]));

% Example image in gray tones and unit-8 image type
Image1 = rgb2gray(Imag1);
Image2 = rgb2gray(Imag2);


% Difraction mask in unit-8 image type
msk = uint8(msk);

mask1 = msk.*Image1;
mask2 = msk.*Image2;


% iamge visualizated trought the LCD simulated
ImageMask1 = mask1;
ImageMask2 = mask2;


figure;
imshow(ImageMask1); title('Image displayed through simulated CCD, Image test 1');

figure;
imshow(ImageMask2); title('Image displayed through simulated CCD, Image test 2');

% Example image in gray tones
figure;
imshow(Image1); title('Image test 1 in gray tones');

figure;
imshow(Image2); title('Image test 1 in gray tones');


% Gray tones profile in the center line of the LCD with image example
cut = mask1( (T/2), : );
if max(cut) == 0
    Cut = mask1( (T/2)+nnn, : );
else 
    Cut = mask1( (T/2), : );
end
x = linspace(1,AcPix,T);
figure;
plot(x,Cut); title('Gray tones distribution at image center');
xlabel('Pixel Number'); ylabel('Gray Tones [0,255]');




%%      Section III
%%%%%%%
%%%%%%% Real Mask Construction
%%%%%%% Mask plane in adition the difractive mask
NN = TT/AcPix;
TTT = int16(fix(T*NN));
Mask0 = zeros(TTT,TTT);
Mask1 = Mask0;
aa = int16(fix(TTT/2));
bb = int16(fix(T/2)); 
%for i = 1:T
Mask1 = uint8(Mask1);
    Mask1( ((aa - bb)):((aa + bb - 1)),((aa - bb)):((aa + bb -1)) ) = ImageMask1;
Mask2 = uint8(Mask0);
    Mask2( ((aa - bb)):((aa + bb - 1)),((aa - bb)):((aa + bb -1)) ) = ImageMask2;

    Mask0( ((aa - bb)):((aa + bb - 1)),((aa - bb)):((aa + bb -1)) ) = DiffMask;


figure;
imshow(Mask1); title('Diffraction Mask, LCD simutaled with Image test 1');

figure;
imshow(Mask2); title('Diffraction Mask, LCD simutaled with Image test 2');


figure;
imshow(Mask0); title('Diffraction Mask without image, only LCD');





%%     Section IV 
%%%%%%
%%%%%% Diffractive patters applied the FFT2



%%%%%%% Diffraction paterns Build, difractive mask with just one square
%%%%%%% aperture. Size aperture PT in milimeters or nnn in matlab pixels
OneMask = zeros(TTT,TTT);

% Wave lenght in millimeters
lambda = 632*(10^(-9));
% K Wave number
K = (2*pi)/lambda;
% Distance from diffractive aperture to obcervation screen, in millimeters 
Z0 = 1000;
R = Z0;
% Fraunhoffer coeficient integral
aaa = (lambda*R);
akz = exp(aaa*1i)/(lambda*R);
Akz = akz*1i;

% Diffractive mas, one square aperture
cc = nnn/2;
    OneMask(  (aa - cc):(aa + cc - 1),(aa - cc):(aa + cc - 1) ) = 1;

% Fraunhoffer integral, FrauU    
FrauU = Akz*fft2( OneMask );
IntFrauU = real( FrauU .* conj( FrauU ) );

%%%%%%%
DiffraFrauU = fftshift( IntFrauU );

figure;
imagesc( DiffraFrauU );
axis( 'equal' ); axis( 'off' );
title( 'Square aperture fraunhoffer diffraction pattern with FFT2' );



%%%%%%% Diffraction paterns Build, difractive mask with just two square
%%%%%%% two squares apertures. Size aperture PT in milimeters or nnn in matlab pixels
TwoMask = zeros(TTT,TTT);

% Diffractive mas, one square aperture
    TwoMask(  (aa - cc):(aa + cc - 1),(aa - cc - nnn):(aa - cc - 1) ) = 1;
    TwoMask(  (aa - cc):(aa + cc - 1),(aa + cc ):(aa + cc + nnn - 1) ) = 1;

% Fraunhoffer integral, FrauU    
FrauTwoU = Akz*fft2( TwoMask );
IntFrauTwoU = real( FrauTwoU .* conj( FrauTwoU ) );

%%%%%%%
DiffraFrauTwoU = fftshift( IntFrauTwoU );

figure;
imagesc( DiffraFrauTwoU );
axis( 'equal' ); axis( 'off' );
title( 'Two simetric Square apertures fraunhoffer diffraction pattern with FFT2' );




%%%%%%% Diffraction paterns Build, difractive mask with just four square
%%%%%%% apertures 2x2. Size aperture PT in milimeters or nnn in matlab pixels
FourMask = zeros(TTT,TTT);

% Diffractive mas, one square aperture
ccc = (3*nnn)/2;


mskfour = zeros(3*nnn, 3*nnn);
% first step, row basis
    mskfour(1:nnn,1:nnn) = 1;
    mskfour(((2*nnn)+1):3*nnn,1:nnn) = 1;
    mskfour(1:nnn,((2*nnn)+1):3*nnn) = 1;
    mskfour(((2*nnn)+1):3*nnn,((2*nnn)+1):3*nnn) = 1;

FourMask( (aa-ccc):((aa+ccc) - 1),(aa-ccc):((aa+ccc) - 1) ) =  mskfour;  

% Fraunhoffer integral, FrauU    
FourFrauU = Akz*fft2( FourMask );
IntFourFrauU = real( FourFrauU .* conj( FourFrauU ) );

%%%%%%%
DiffraFourFrauU = fftshift( IntFourFrauU );


figure;
imagesc( DiffraFourFrauU );
axis( 'equal' ); axis( 'off' );
title( 'Four 2x2 Square aperture fraunhoffer diffraction pattern with FFT3' );





%%%%%%% Diffraction paterns Build, difractive mask with just nine square
%%%%%%% apertures 3x3. Size aperture PT in milimeters or nnn in matlab pixels
NineMask = zeros(TTT,TTT);


% Diffractive mas, one square aperture
cccc = (5*nnn)/2;


msknine = zeros(5*nnn, 5*nnn);
% first step, row basis
    msknine(1:nnn,1:nnn) = 1;
    msknine(1:nnn,((2*nnn)+1):3*nnn) = 1;
    msknine(1:nnn,((4*nnn)+1):5*nnn) = 1;
    
    msknine(((2*nnn)+1):3*nnn,1:nnn) = 1;
    msknine(((2*nnn)+1):3*nnn,((2*nnn)+1):3*nnn) = 1;
    msknine(((2*nnn)+1):3*nnn,((4*nnn)+1):5*nnn) = 1;
    
    msknine(((4*nnn)+1):5*nnn,1:nnn) = 1;
    msknine(((4*nnn)+1):5*nnn,((2*nnn)+1):3*nnn) = 1;
    msknine(((4*nnn)+1):5*nnn,((4*nnn)+1):5*nnn) = 1;
    

NineMask( (aa-cccc):((aa+cccc) - 1),(aa-cccc):((aa+cccc) - 1) ) =  msknine;  

% Fraunhoffer integral, FrauU    
NineFrauU = Akz*fft2( NineMask );
IntNineFrauU = real( NineFrauU .* conj( NineFrauU ) );

%%%%%%%
DiffraNineFrauU = fftshift( IntNineFrauU );

figure;
imagesc( DiffraNineFrauU );
axis( 'equal' ); axis( 'off' );
title( 'Nine 3x3 Square aperture fraunhoffer diffraction pattern with FFT2' );




%% 
%%%%% Seccion V
%%%%% Diffractive patterns following the calculations in the document as
%%%%% "supplementary material"
%%%%% All calculations are in centimeters
% Aberture
a = 6*(10^(-4));
b = 6*(10^(-4));
% z distance, distance between the diffraction window and the observation plane
L = 10;
c = 512;
% LCSLM (LCD) Pixel number N




N = 1;
% LCD size D = 2d
d = 2.54/2;
y = linspace(-d,d,c);
x=imresize(y,[1 c]);
[X,Y] = meshgrid(x,x);
% wavelength
LO = 632*(10^(-7));
z = 0.1;
aa = ( a*(LO) ) / ( 2*(pi)*L );
Aa = (aa.*aa);
sincx = sinc(X*(  a*((pi)/(L*LO)) ));
Sincx = sincx.*sincx;
sincy = sinc(Y*(  a*((pi)/(L*LO)) ));
Sincy = sincy.*sincy;
I1x1 = Aa.*Sincx.*Sincy;

figure;
mesh(X, Y, I1x1), colormap jet
xlabel('x [cm]'); ylabel('y [cm]'); zlabel('I(x,y) [W/m2]');
title('LCSLM (LCD) 1x1 pixels Diffractive patern');





% considering the N-th exponent of the geometric series on Sab(x) or Sab(y)

bbx = X.*(( (2*(pi))/(z*LO) )*(a+b));
BBx = complex(0,bbx);
bby = Y.*(( (2*(pi))/(z*LO) )*(a+b));
BBy = complex(0,bby);
% considering the N-th exponent of the geometric series on Sab(x) or Sab(y)
N = 2;
n= N - 1;

Mdiffx = zeros(c,c,n);
BMatx = zeros(c,c,n);
Mdiffy = zeros(c,c,n);
BMaty = zeros(c,c,n);

for j=1:N
    BMatx(:,:,j) = (j-1).*BBx;
    Mdiffx(:,:,j) = exp(-(  BMatx(:,:,j)  ));
    BMaty(:,:,j) = (j-1).*BBy;
    Mdiffy(:,:,j) = exp(-(  BMaty(:,:,j)  ));
end

sabXn1 = sum(Mdiffx,3);
SabXn1 = real( sabXn1 .* conj( sabXn1 ) );
sabYn1 = sum(Mdiffy,3);
SabYn1 = real( sabYn1 .* conj( sabYn1 ) );

% 2x2 pixel simulated LCD intensity
I2x2 = Aa.*Sincx.*Sincy.*SabXn1.*SabYn1;

figure;
mesh(X, Y, I2x2), colormap jet
xlabel('x [cm]'); ylabel('y [cm]'); zlabel('I(x,y) [W/m2]');
title('LCSLM (LCD) 2x2 pixels Diffractive patern');





% considering the N-th exponent of the geometric series on Sab(x) or Sab(y)
N = 3;
n= N - 1;

Mdiffx = zeros(c,c,n);
BMatx = zeros(c,c,n);
Mdiffy = zeros(c,c,n);
BMaty = zeros(c,c,n);

for j=1:N
    BMatx(:,:,j) = (j-1).*BBx;
    Mdiffx(:,:,j) = exp(-(  BMatx(:,:,j)  ));
    BMaty(:,:,j) = (j-1).*BBy;
    Mdiffy(:,:,j) = exp(-(  BMaty(:,:,j)  ));
end

sabXn1 = sum(Mdiffx,3);
SabXn1 = real( sabXn1 .* conj( sabXn1 ) );
sabYn1 = sum(Mdiffy,3);
SabYn1 = real( sabYn1 .* conj( sabYn1 ) );

% 2x2 pixel simulated LCD intensity
I3x3 = Aa.*Sincx.*Sincy.*SabXn1.*SabYn1;

figure;
mesh(X, Y, I3x3), colormap jet
xlabel('x [cm]'); ylabel('y [cm]'); zlabel('I(x,y) [W/m2]');
title('LCSLM (LCD) 3x3 pixels Diffractive patern');





% considering the N-th exponent of the geometric series on Sab(x) or Sab(y)
N = 4;
n= N - 1;

Mdiffx = zeros(c,c,n);
BMatx = zeros(c,c,n);
Mdiffy = zeros(c,c,n);
BMaty = zeros(c,c,n);

for j=1:N
    BMatx(:,:,j) = (j-1).*BBx;
    Mdiffx(:,:,j) = exp(-(  BMatx(:,:,j)  ));
    BMaty(:,:,j) = (j-1).*BBy;
    Mdiffy(:,:,j) = exp(-(  BMaty(:,:,j)  ));
end

sabXn1 = sum(Mdiffx,3);
SabXn1 = real( sabXn1 .* conj( sabXn1 ) );
sabYn1 = sum(Mdiffy,3);
SabYn1 = real( sabYn1 .* conj( sabYn1 ) );

% 2x2 pixel simulated LCD intensity
I4x4 = Aa.*Sincx.*Sincy.*SabXn1.*SabYn1;

figure;
mesh(X, Y, I4x4), colormap jet
xlabel('x [cm]'); ylabel('y [cm]'); zlabel('I(x,y) [W/m2]');
title('LCSLM (LCD) 4x4 pixels Diffractive patern');





% considering the N-th exponent of the geometric series on Sab(x) or Sab(y)
%%%%  This is a general case for all N
N = 512;
n= N - 1;

Mdiffx = zeros(c,c,n);
BMatx = zeros(c,c,n);
Mdiffy = zeros(c,c,n);
BMaty = zeros(c,c,n);

for j=1:N
    BMatx(:,:,j) = (j-1).*BBx;
    Mdiffx(:,:,j) = exp(-(  BMatx(:,:,j)  ));
    BMaty(:,:,j) = (j-1).*BBy;
    Mdiffy(:,:,j) = exp(-(  BMaty(:,:,j)  ));
end

sabXn1 = sum(Mdiffx,3);
SabXn1 = real( sabXn1 .* conj( sabXn1 ) );
sabYn1 = sum(Mdiffy,3);
SabYn1 = real( sabYn1 .* conj( sabYn1 ) );

% 2x2 pixel simulated LCD intensity
I512x512 = Aa.*Sincx.*Sincy.*SabXn1.*SabYn1;

figure;
mesh(X, Y, I512x512), colormap jet
xlabel('x [cm]'); ylabel('y [cm]'); zlabel('I(x,y) [W/m2]');
title('LCSLM (LCD) 512x512 pixels Diffractive patern');







%%   Section VI
%%%% Super Gaussian Profiles and Super-Gaussian fringes ruling

c = 1024;
y = -(c/2):1/c:(c/2);
x=imresize(y,[1 c]);
[X,Y] = meshgrid(x,x);
R0 = 1;
r = (sqrt( (X.*X) + (Y.*Y) ));
a = 1/12;            % Valor solo para el ejemplo de 50 franjas SG
                     %este termino determina el ancho de las franjas, su 
                     %ancho debe ser el valor t, pero al incrementar a se 
                     %disminulle el ancho y al aumentar a se incrementa el
                     %ancho de las ffranjas SG
                     
sigma = 1/(sqrt(a));
Wb =  2*(sigma^2);

gamma = 100;

SG = R0*(exp( - ((r.^2)./(Wb) ).^gamma ) );

figure;
imshow(SG)
title('Super-Gaussian profile');



SGx = R0*(exp( - ((X.^2)./(Wb) ).^gamma ) );

figure;
imshow(SGx)
title('Super-Gaussian fringe, pattern basis');


% Super Gaussian center in a X0
% center at
t = 10.24;       % solo para el ejemplo de 50 franjas

%t = c/(2*Num);   %expresion general que determina el ancho de las franjas
                  %segun se incremente o dismiya el ancho de las franjas
                  %respecto del valor de "a"

% Because "c = 1024", but if "c = 512" then it must be satisfy that" t = 5.12"

X0 = (t)-(c/2);

SGx0 = R0*(exp( -  ((( (X-X0).^2 )./(Wb)).^gamma)     ));

figure;
imshow(SGx0)
title('Super-Gaussian fringe, pattern basis at the beginning ');



% Super-Gaussian Periodic patern, as vertica Ronchi rruling

Num = 50;
% Number of the Super-Gaussian fringes

SGXk = zeros(c,c,Num);

SGk(:,:,1) = SGx0;

n = Num -1;
for i=1:n
    
    Y0 = t+(i*(2*t)) - (c/2);
    SGk(:,:,i+1) = R0*(exp( -  ((( (X-Y0).^2 )./(Wb)).^gamma)     ));
    
end

% Super-Gaussian rulling
A = sum(SGk,3);

figure;
imshow(A)
title('Super-Gaussian fringe ruling pattern as Ronchi ruling');


% profile plot
Cut = A(512,:);

figure
plot(X(512,:), Cut)
title('Super-Gaussian fringe ruling profile');



% Difraction mask in unit-8 image type with Super-Gaussian profile

[dd ee] = size(msk);

AA = imresize(A, [dd ee]);

MSK = double(msk);
mask4 = MSK.*AA;


% Imamge visualizated trought the LCD simulated
ImageMask4 = mask4;


figure;
imshow(ImageMask4); title('Image displayed through simulated CCD, Image test 1');



%%%%%%% Mask plane in adition the difractive mask
 
Mask3 = zeros(TTT,TTT);

    Mask3( ((aa - bb)):((aa + bb - 1)),((aa - bb)):((aa + bb -1)) ) = ImageMask4;


figure;
imshow(Mask3); title('Diffraction Mask, LCD simutaled with Image test 1');
