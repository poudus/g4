<?php
$sid = $_GET["sid"];

	$connstring = 'dbname=g4 user=g4';
	$conn = pg_connect($connstring) or die ("could not connect");
	$query = 'select * from g4s.frame where sid = ' . $sid . ' order by epoch, pid asc';
	$rs = pg_query($conn, $query);

	while ($row = pg_fetch_row($rs)) {
		echo json_encode($row);
	}
	pg_close($conn);
?>
