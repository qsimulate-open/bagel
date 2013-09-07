<?php
$authors = "https://raw.github.com/nubakery/bagel/master/doc/authors.json";
$file = file_get_contents($authors);
$json = json_decode($file, true);
echo '<ul>';
foreach ($json as $name) {
    echo '<li style="margin-top:10px"><strong>'.$name['name'].' '.$name['surname'] .'</strong> ('.$name['start'].'&ndash;'.$name['end'].') ';
    echo ($name['affiliation']==""?'&nbsp;':$name['affiliation']).'<br />';
    echo ($name['email']==""?'&nbsp;':$name['email']).'</li>';
}
echo '</ul>';
echo '<div style="margin-top:30px">Generated from bagel/doc/<a href="'.$authors.'">authors.json</a></div>';
?>
