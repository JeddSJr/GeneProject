import { Component, ViewChild, ViewContainerRef} from '@angular/core';

import { AuthService } from '@auth0/auth0-angular';
@Component({
  selector: 'app-home',
  templateUrl: './home.component.html',
  styleUrl: './home.component.css'
})
export class HomeComponent {
  loggedIn = false;


  @ViewChild("viewContainerRef", {read: ViewContainerRef})
    // @ts-ignore
  VCR: ViewContainerRef;

  constructor(public auth: AuthService){
    this.auth.isAuthenticated$.subscribe(authenticated =>{
      this.checkLogin(authenticated);
    })
  }

  checkLogin(authenticated: boolean){
    this.loggedIn = authenticated
  }
}
