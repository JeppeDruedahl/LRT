function RGB = hex2rgb(HEX)

	if ~iscellstr(HEX)
		HEX = cellstr(HEX);
	end
	HEX = upper(HEX);
	HEX = double(char(HEX));
	tf = HEX > 64;
	HEX(~tf) = HEX(~tf) - '0';
	HEX(tf) = HEX(tf) - 'A' + 10;

	HEX = reshape(HEX.',2,[]);
	HEX(1,:) = 16*HEX(1,:);
	RGB = sum(HEX,1);
	RGB  = reshape(RGB,3,[]);
	RGB = RGB.';

end