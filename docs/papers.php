<?php
$papers = "https://raw.github.com/nubakery/bagel/master/doc/papers.json";
$file = file_get_contents($papers);
$json = json_decode($file, true);
$i = 0;
foreach ($json as $year) {
   echo '<div>'.$year['year'].'</div>';
   echo '<div><ul>';
   foreach ($year['references'] as $pub) {
     $i = $i+1;
     echo '<li style="margin-top:10px">';
     echo $pub['author'].',<br />';
     echo '<i>'.$pub['journal'].'</i> <b>'.$pub['volume'].'</b>, '.$pub['pages'].' ('.$year['year'].'),<br />';
     echo '<a href="http://dx.doi.org/'.$pub['doi'].'">&ldquo;'.$pub['title'].'&rdquo;</a>';
     echo '</li>';
   }
   echo '</ul></div>';
}
echo '<div style="margin-top:30px">Generated from bagel/doc/<a href="'.$papers.'">papers.json</a></div>';
?>
