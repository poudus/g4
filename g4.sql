
create schema g4s;

create sequence g4s.simulation_sequo start 1;

create table g4s.simulation (
id int primary key,
ts bigint not null,
epochs int not null,
frames int not null,
gravity numeric(6,2) not null,
rebound numeric(6,2) not null,
transfer numeric(6,2) not null,
mmin int not null,
mmax int not null,
vmin int not null,
vmax int not null,
bc int not null,
pc int not null,
views int not null,
likes int not null
);

create table g4s.frame (
sid int not null,
epoch int not null,
pid int not null,
x int not null,
y int not null,
radius numeric(6,2) not null,
color varchar(8) not null
);

