package main

import (
	"github.com/fsnotify/fsnotify"
	"log"
	"os"
	"time"
)

// IsDirectory - Returns a true/false if file is a directory
func IsDirectory(path string) bool {
	fileInfo, err := os.Stat(path)
	if err != nil {
		log.Fatal(err)
	}
	return fileInfo.IsDir()
}

func main() {
	watcher, err := fsnotify.NewWatcher()
	if err != nil {
		log.Fatal(err)
	}
	defer watcher.Close()
	done := make(chan bool)

	go func() {
		for {
			select {
			case event, ok := <-watcher.Events:
				if !ok {
					return
				}
				log.Println("event:", event)
				if IsDirectory(event.Name) == true {
					log.Println("Adding new directory to watch: ", event.Name)
					watcher.Add(event.Name)
				}
				if event.Op&fsnotify.Write == fsnotify.Write {
					log.Println("modified file:", event.Name)
				}
			case err, ok := <-watcher.Errors:
				if !ok {
					return
				}
				log.Println("error:", err)
			case <-time.After(10 * time.Second):
				log.Println("TIMEOUT!")
				close(done)
			}
		}
	}()

	err = watcher.Add(os.Args[1])
	if err != nil {
		log.Fatal(err)
	}
	<-done
}
