import { Component } from '@angular/core';
import { AuthService } from '@auth0/auth0-angular';
import { DOCUMENT } from '@angular/common';
import { Inject} from '@angular/core';
@Component({
  selector: 'app-auth-button',
  templateUrl: './auth-button.component.html',
  styleUrl: './auth-button.component.css'
})
export class AuthButtonComponent {
  constructor(@Inject(DOCUMENT) public document: Document, public auth: AuthService) {}
  
  loginWithRedirect(){
    this.auth.loginWithRedirect()
    console.log(this.auth.isAuthenticated$)
  }
  logout(){
    this.auth.logout({ logoutParams: { returnTo: document.location.origin } })
   }
}
