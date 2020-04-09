function reg_img = register_img(raw_img)

sz=size(raw_img);
frtoavg=20;
Reg_stack = zeros(sz(1),sz(2),1,'uint16');
reg_img = zeros(sz(1),sz(2),frtoavg,'uint16');
xyshift = zeros(frtoavg, 2);
avg=raw_img(:,:,1);
for fr=1:sz(3);
    [output, Greg ] = dftregistration(fft2(avg(50:end-50,50:end-50)),fft2(raw_img(50:end-50,50:end-50,fr)),1);%IG cropping edges
    xyshift(fr, :)  = [output(3) output(4)];
    reg_error(fr)=output(1);
    transvec=[output(3),output(4)];%MFI
    se=translate(strel(1),transvec);%MFI
    Reg_stack(:,:,1) = uint16(imdilate(raw_img(:,:,fr),se));%MFI
    reg_img(:,:,fr) = Reg_stack(:,:,1);
    if fr<= frtoavg
        avg=mean(reg_img,3);
    end
end