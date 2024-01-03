function [outputArg1] = dec2sex(dec)
%Decimal 2 sexagesimal degrees
deg = fix(dec);
primes = fix((dec - deg)*60);
seconds = ((dec - deg)*60  - primes)*60;

outputArg1 = [deg; primes; seconds];

end